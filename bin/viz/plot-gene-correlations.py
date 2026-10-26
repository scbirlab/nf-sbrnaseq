
#!/usr/bin/env python
"""Gene correlations and functional annotations: COG, EC, KEGG."""

import argparse
from collections import Counter
from itertools import cycle
from pathlib import Path

import anndata as ad
from carabiner import colorblind_palette, print_err
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.sparse import issparse
from scipy.spatial.distance import squareform
import umap

from viz_utils import heatmap


ANNOTATIONS = {
    "cog": "COG_category",
    "ec": "EC",
    "kegg": "KEGG_Pathway",
}


def memberships(value, kind):
    """Parse functional memberships without losing multi-annotations."""
    if pd.isna(value) or str(value).strip() in ("", "-", "."):
        return ()

    value = str(value).strip()
    terms = list(value) if kind == "cog" else value.split(",")

    return tuple(sorted({
        term.strip() for term in terms
        if term.strip() and term.strip() not in ("-", ".")
    }))


def primary_groups(values, kind):
    """Assign one group per gene using global group prevalence."""
    parsed = [memberships(value, kind) for value in values]
    counts = Counter(
        term for terms in parsed for term in terms
    )

    groups = [
        min(terms, key=lambda term: (-counts[term], term))
        if terms else "Unclassified"
        for terms in parsed
    ]

    return parsed, groups



def read_annotations(path, genes):
    """Join functional annotations using AnnData composite gene IDs."""
    annotations = pd.read_csv(
        path, sep="\t", dtype=str, keep_default_na=False
    )

    gene_columns = [
        "genome_accession", "assembly", "chr", "gene_id"
    ]
    required = set(gene_columns) | set(ANNOTATIONS.values())

    missing = required - set(annotations.columns)
    if missing:
        raise ValueError(
            f"Missing annotation columns: {sorted(missing)}"
        )

    annotations["composite_id"] = annotations[
        gene_columns
    ].agg(":".join, axis=1)

    if annotations["composite_id"].duplicated().any():
        raise ValueError("Duplicate composite gene identifiers")

    lookup = annotations.set_index("composite_id")
    matched = pd.Index(genes).isin(lookup.index)
    n_matched = int(matched.sum())

    print_err(
        f"Annotation join: {n_matched}/{len(genes)} genes "
        f"({n_matched / len(genes):.1%})"
    )

    if n_matched == 0:
        raise ValueError(
            "No gene identifiers matched functional annotations"
        )

    result = lookup.reindex(genes)[
        list(ANNOTATIONS.values())
    ].fillna("")

    result.index.name = "gene_id"
    return result


def expression_matrix(adata, min_cells, max_genes=None):
    """Select variable genes using detection in the raw count layer."""
    if "counts" not in adata.layers:
        raise ValueError("Missing counts layer")

    detected = np.asarray(
        (adata.layers["counts"] > 0).sum(axis=0)
    ).ravel()

    eligible = np.flatnonzero(detected >= min_cells)
    if len(eligible) < 2:
        raise ValueError("Fewer than two eligible genes")

    x = adata.X[:, eligible]
    x = x.toarray() if issparse(x) else np.asarray(x)

    chosen = np.argsort(
        -x.var(axis=0), 
        kind="stable",
    )
    
    if max_genes is not None and isinstance(max_genes, int):
        chosen = chosen[:max_genes]

    frame = pd.DataFrame(
        x[:, chosen],
        index=adata.obs_names,
        columns=adata.var_names[eligible[chosen]],
    )

    frame = frame.loc[:, frame.nunique() > 1]
    if frame.shape[1] < 2:
        raise ValueError("Fewer than two variable genes")

    return frame


def correlations(frame):
    """Compute correlations, distances and clustering once."""
    corr = frame.corr(method="spearman")
    values = corr.to_numpy()

    if not np.isfinite(values).all():
        raise ValueError("Non-finite gene correlations")

    distances = np.clip(1 - values, 0, 2)
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0)

    hierarchy = (
        linkage(
            squareform(distances, checks=True),
            method="average",
        )
        if len(corr) > 2 else None
    )

    return corr, distances, hierarchy


def gene_embedding(distances):
    """One deterministic embedding shared by all classifications."""
    n = distances.shape[0]

    if n < 4:
        return np.full((n, 2), np.nan)

    return umap.UMAP(
        # n_neighbors=min(15, n - 1),
        metric="euclidean",
        min_dist=.1,
        random_state=42,
    ).fit_transform(distances)


def group_colors(groups, top):
    """Colour all groups, prioritising those shown on UMAP."""
    names = sorted(set(groups) - {"Unclassified"})

    ordered = list(top) + [
        name for name in names if name not in top
    ]

    palette = cycle(colorblind_palette())
    colors = {
        name: next(palette) for name in ordered
    }
    colors["Unclassified"] = "lightgrey"

    return colors


def plot_gene_umap(
    coords, 
    groups, 
    colors, 
    top,
    kind, 
    filename, 
    fmt
):
    """Show eight most populated functional groups."""
    if not np.isfinite(coords).all():
        print_err(f"Skipping {kind} UMAP: too few genes")
        return

    fig, ax = grid(figsize=(8, 6))
    groups = np.asarray(groups)

    background = ~np.isin(groups, top)
    ax.scatter(
        *coords[background].T,
        color="lightgrey",
        s=12,
        alpha=0.65,
        rasterized=True,
        label="Other / unclassified",
    )

    for group in top:
        ma = groups == group
        ax.scatter(
            *coords[mask].T,
            color=colors[group],
            s=20,
            rasterized=True,
            label=f"{group} (n={mask.sum()})",
        )

    ax.set(
        xlabel="Gene UMAP 1",
        ylabel="Gene UMAP 2",
        title=f"Gene correlations — {kind.upper()}",
    )
    ax.legend(
        title=ANNOTATIONS[kind],
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=8,
    )

    fig.savefig(
        f"{filename}.{kind}.gene-umap.{fmt}",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clustered", required=True)
    parser.add_argument("--annotations", default=None)
    parser.add_argument("--output", default="correlations")
    parser.add_argument(
        "--min-cells-per-gene", type=int, default=3
    )
    parser.add_argument("--max-genes", type=int, default=None)
    parser.add_argument("--top-groups", type=int, default=8)
    parser.add_argument("--plot-format", default="png")
    args = parser.parse_args()

    if args.max_genes is not None:
        if args.max_genes < 2:
            parser.error("--max-genes must be >= 2")
        if args.top_groups < 1:
            parser.error("--top-groups must be >= 1")

    prefix = Path(args.output)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    rna = ad.read_h5ad(args.clustered)
    frame = expression_matrix(
        rna, 
        args.min_cells_per_gene, 
        args.max_genes,
    )

    corr, distances, hierarchy = correlations(frame)
    corr.to_csv(f"{prefix}.matrix.csv")

    coords = gene_embedding(distances)
    df = pd.DataFrame(
        coords,
        index=corr.index,
        columns=["UMAP_1", "UMAP_2"],
    )

    size = max(7, min(28, len(corr) * .15))

    # Without functional annotations, produce one plain
    # correlation heatmap and one uncoloured gene UMAP.
    if not args.annotations:
        print_err("No eggNOG annotations: plotting plain correlations")

        hm = heatmap(
            corr,
            color_label="Spearman $\\rho$",
            filename=f"{prefix}.heatmap",
            output_format=args.plot_format,
            row_linkage=hierarchy,
            col_linkage=hierarchy,
            cmap="vlag",
            vmin=-1,
            vmax=1,
            center=0,
            figsize=(size, size),
            mask_missing=False,
        )

        if np.isfinite(coords).all():
            fig, ax = grid(aspect_ratio=1.2)
            ax.scatter(
                *coords.T,
                s=10.,
                rasterized=True,
            )
            ax.set(
                xlabel="Gene UMAP 1",
                ylabel="Gene UMAP 2",
                title="Gene expression correlations",
            )
            _figsaver(
                fig=fig,
                name=f"{prefix}.gene-umap",
                df=df,
            )
        else:
            print_err("Skipping gene UMAP: fewer than four genes")

        print_err(f"Analysed {len(corr)} genes")
        return None

    # Functional annotations: COG, EC and KEGG.
    metadata = read_annotations(
        args.annotations, 
        corr.index,
    )
    metadata["UMAP_1"], metadata["UMAP_2"] = coords[:, :2].T

    for kind, field in ANNOTATIONS.items():
        parsed, groups = primary_groups(
            metadata[field], 
            kind,
        )

        counts = Counter(
            g for g in groups if g != "Unclassified"
        )
        top = [
            group for group, _ in sorted(
                counts.items(),
                key=lambda x: (-x[1], x[0]),
            )[:args.top_groups]
        ]

        colors = group_colors(groups)

        metadata[f"{kind}_memberships"] = [
            ",".join(terms) for terms in parsed
        ]
        metadata[f"{kind}_group"] = groups
        metadata[f"{kind}_heatmap_color"] = [
            colors[group] for group in groups
        ]
        metadata[f"{kind}_umap_group"] = [
            group if group in top else "Other / unclassified"
            for group in groups
        ]
        metadata[f"{kind}_umap_color"] = [
            colors[group] if group in top else "lightgrey"
            for group in groups
        ]

        strips = pd.Series(
            metadata[f"{kind}_heatmap_color"],
            index=corr.index,
            name=kind.upper(),
        )

        hm = heatmap(
            corr,
            color_label="Spearman correlation",
            filename=f"{prefix}.{kind}.heatmap",
            output_format=args.plot_format,
            row_colors=strips,
            col_colors=strips,
            row_linkage=hierarchy,
            col_linkage=hierarchy,
            cmap="vlag",
            vmin=-1,
            vmax=1,
            center=0,
            figsize=(size, size),
            mask_missing=False,
        )

        plot_gene_umap(
            coords, groups, colors, top,
            kind, prefix, args.plot_format,
        )

        pd.DataFrame([
            {
                "group": group,
                "color": color,
                "highlighted_on_umap": group in top,
                "n_genes": counts.get(group, 0),
            }
            for group, color in sorted(colors.items())
        ]).to_csv(
            f"{prefix}.{kind}.colors.csv",
            index=False,
        )

        print_err(
            f"{kind.upper()}: {len(counts)} groups; "
            f"{len(top)} highlighted on UMAP"
        )

    metadata.to_csv(f"{prefix}.annotations.csv")
    print_err(f"Analysed {len(corr)} genes")


if __name__ == "__main__":
    main()
