
#!/usr/bin/env python3
"""Gene correlations with COG, EC and KEGG annotation heatmaps."""

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
import seaborn as sns


ANNOTATIONS = {
    "cog": "COG_category",
    "ec": "EC",
    "kegg": "KEGG_Pathway",
}


def memberships(value, kind):
    """Parse all memberships, preserving multiple annotations."""
    if pd.isna(value) or str(value).strip() in ("", "-", "."):
        return ()

    value = str(value).strip()
    terms = list(value) if kind == "cog" else value.split(",")
    return tuple(sorted({
        term.strip() for term in terms if term.strip()
    }))


def assign_groups(values, kind):
    """Choose one primary group per gene deterministically.

    Multi-annotated genes are assigned to their most prevalent
    membership among the plotted genes. Lexical order breaks ties.
    """
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


def annotation_table(path, genes):
    """Match composite AnnData gene IDs to GFF-derived annotations."""
    annotations = pd.read_csv(
        path, sep="\t", dtype=str, keep_default_na=False
    )

    required = {
        "contig", "gene_id", *ANNOTATIONS.values()
    }
    missing = required - set(annotations.columns)
    if missing:
        raise ValueError(f"Missing annotation columns: {missing}")

    if annotations.duplicated(["contig", "gene_id"]).any():
        raise ValueError("Duplicate (contig, gene_id) annotation keys")

    lookup = annotations.set_index(["contig", "gene_id"])

    keys = [str(gene).split(":", 3) for gene in genes]
    if any(len(key) != 4 for key in keys):
        raise ValueError(
            "Expected genome_accession:assembly:chr:gene_id"
        )

    rows = []
    matched = 0

    for genome, assembly, contig, gene_id in keys:
        key = (contig, gene_id)
        if key in lookup.index:
            record = lookup.loc[key]
            matched += 1
            rows.append({
                field: record[field]
                for field in ANNOTATIONS.values()
            })
        else:
            rows.append({
                field: "" for field in ANNOTATIONS.values()
            })

    print_err(
        f"Annotation join: {matched}/{len(genes)} genes "
        f"({matched / len(genes):.1%})"
    )

    result = pd.DataFrame(rows, index=genes)
    result.index.name = "gene_id"
    return result


def expression_matrix(adata, min_cells, max_genes):
    """Select variable, sufficiently detected genes."""
    if "counts" not in adata.layers:
        raise ValueError("Missing raw counts layer")

    detected = np.asarray(
        (adata.layers["counts"] > 0).sum(axis=0)
    ).ravel()

    eligible = np.flatnonzero(detected >= min_cells)
    if len(eligible) < 2:
        raise ValueError("Fewer than two qualifying genes")

    expression = adata.X[:, eligible]
    if issparse(expression):
        expression = expression.toarray()
    else:
        expression = np.asarray(expression)

    variance = expression.var(axis=0)
    chosen = np.argsort(
        -variance, kind="stable"
    )[:max_genes]

    frame = pd.DataFrame(
        expression[:, chosen],
        index=adata.obs_names,
        columns=adata.var_names[eligible[chosen]],
    )

    frame = frame.loc[:, frame.nunique() > 1]

    if frame.shape[1] < 2:
        raise ValueError("Fewer than two variable genes")

    return frame


def correlation_hierarchy(correlation):
    """Return one shared average-linkage gene hierarchy."""
    distances = np.clip(
        1 - correlation.to_numpy(), 0, 2
    )
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0)

    if len(correlation) <= 2:
        return None

    return linkage(
        squareform(distances, checks=True),
        method="average",
    )


def plot_heatmap(correlation, hierarchy, groups,
                 colors, filename, kind, fmt):
    """Draw correlation matrix with matching row/column strips."""
    row_colors = pd.Series(
        [colors[group] for group in groups],
        index=correlation.index,
        name=kind.upper(),
    )

    size = max(7, min(28, len(correlation) * 0.15))

    hm = sns.clustermap(
        correlation,
        row_linkage=hierarchy,
        col_linkage=hierarchy,
        row_cluster=hierarchy is not None,
        col_cluster=hierarchy is not None,
        row_colors=row_colors,
        col_colors=row_colors,
        cmap="vlag",
        vmin=-1,
        vmax=1,
        center=0,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Spearman correlation"},
        figsize=(size, size),
    )

    hm.ax_heatmap.tick_params(
        axis="both", labelsize=5
    )

    hm.savefig(
        f"{filename}.{kind}.heatmap.{fmt}",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(hm.fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clustered", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--output", default="correlations")
    parser.add_argument(
        "--min-cells-per-gene", type=int, default=3
    )
    parser.add_argument(
        "--max-genes", type=int, default=150
    )
    parser.add_argument("--plot-format", default="png")
    args = parser.parse_args()

    if args.max_genes < 2:
        parser.error("--max-genes must be >= 2")

    prefix = Path(args.output)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.clustered)
    frame = expression_matrix(
        adata, args.min_cells_per_gene, args.max_genes
    )

    corr = frame.corr(method="spearman")
    if not np.isfinite(corr.to_numpy()).all():
        raise ValueError("Non-finite gene correlations")

    corr.to_csv(f"{prefix}.matrix.csv")
    hierarchy = correlation_hierarchy(corr)

    annotation = annotation_table(
        args.annotations, corr.index
    )

    # Same gene order and hierarchy in every heatmap.
    for kind, column in ANNOTATIONS.items():
        parsed, groups = assign_groups(
            annotation[column], kind
        )

        unique_groups = sorted(
            set(groups) - {"Unclassified"}
        )

        palette = cycle(colorblind_palette())
        colors = {
            group: next(palette)
            for group in unique_groups
        }
        colors["Unclassified"] = "lightgrey"

        annotation[f"{kind}_memberships"] = [
            ",".join(terms) for terms in parsed
        ]
        annotation[f"{kind}_group"] = groups
        annotation[f"{kind}_color"] = [
            colors[group] for group in groups
        ]

        plot_heatmap(
            corr, hierarchy, groups, colors,
            prefix, kind, args.plot_format,
        )

        pd.DataFrame([
            {"group": group, "color": color}
            for group, color in sorted(colors.items())
        ]).to_csv(
            f"{prefix}.{kind}.colors.csv",
            index=False,
        )

        print_err(
            f"{kind.upper()}: {len(unique_groups)} groups, "
            f"{sum(g != 'Unclassified' for g in groups)} "
            "annotated genes"
        )

    annotation.to_csv(f"{prefix}.annotations.csv")
    print_err(
        f"Generated three heatmaps for {len(corr)} genes"
    )


if __name__ == "__main__":
    main()
