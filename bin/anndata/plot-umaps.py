
#!/usr/bin/env python
"""Plot RNA and multimodal UMAP feature grids."""

import argparse
from pathlib import Path

import anndata as ad
from carabiner.mpl import add_legend, grid, figsaver
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import issparse


def dense_column(matrix, index):
    value = matrix[:, index]
    if issparse(value):
        return value.toarray().ravel()
    return np.asarray(value).ravel()


def row_sums(matrix):
    return np.asarray(matrix.sum(axis=1)).ravel()


def biotype_groups(adata):
    if "gene_biotype" not in adata.var:
        raise ValueError("Missing gene_biotype annotations")

    biotype = (
        adata.var["gene_biotype"]
        .astype(str)
        .fillna("")
        .str.lower()
    )

    trna = biotype.eq("trna")
    tmrna = biotype.isin(["tmrna", "transfer-messenger rna"])
    rrna = biotype.eq("rrna")

    other_ncrna = (
        biotype.str.contains("ncrna|non_coding|non-coding|srna|antisense")
        | biotype.isin(["rna", "scrna", "snrna"])
    ) & ~(trna | tmrna | rrna)

    return {
        "trna_tmrna": adata.var_names[trna | tmrna].tolist(),
        "ncrna": adata.var_names[other_ncrna].tolist(),
    }


def feature_label(adata, feature):
    row = adata.var.loc[feature]
    label = str(row.get("Name", ""))
    if label and label.lower() != "nan" and label != feature:
        return f"{label}\t{feature}"
    return str(feature)


def plot_grid(coords, frame, columns, filename, fmt, ncol=6):
    if not columns:
        return

    saver = figsaver(format=fmt)

    # Paginate so large tRNA catalogues remain readable.
    per_page = 24

    for offset in range(0, len(columns), per_page):
        page = columns[offset:offset + per_page]
        page_number = offset // per_page + 1

        fig, axes = grid(
            ncol=min(ncol, len(page)),
            nrow=int(np.ceil(len(page) / ncol)),
            aspect_ratio=1.1,
            panel_size=5.,
        )
        axes = np.asarray(axes).ravel()

        for ax, column in zip(axes, page):
            values = frame[column]
            valid = values.notna().to_numpy()

            ax.scatter(
                coords[~valid, 0],
                coords[~valid, 1],
                c="lightgrey",
                s=1.,
                rasterized=True,
            )

            if isinstance(values.dtype, pd.CategoricalDtype) or (
                pd.api.types.is_object_dtype(values.dtype)
                or pd.api.types.is_string_dtype(values.dtype)
            ):
                categories = values[valid].astype(str)
                for category in sorted(categories.unique()):
                    selected = (
                        (values.astype(str) == category).to_numpy()
                        & valid
                    )
                    ax.scatter(
                        coords[selected, 0],
                        coords[selected, 1],
                        s=1.5,
                        label=category,
                        rasterized=True,
                    )
                add_legend(ax)
            else:
                numeric = pd.to_numeric(values, errors="coerce")
                scatter = ax.scatter(
                    coords[valid, 0],
                    coords[valid, 1],
                    c=numeric.to_numpy()[valid],
                    cmap="cividis",
                    s=1.5,
                    rasterized=True,
                )
                fig.colorbar(scatter, ax=ax)

            ax.set(
                title=str(column),
                xlabel="UMAP 1",
                ylabel="UMAP 2",
            )
            ax.set_aspect("equal")

        for ax in axes[len(page):]:
            ax.axis("off")

        saver(
            fig=fig,
            name=f"{filename}.page{page_number:02d}",
            df=frame[page].assign(
                UMAP_1=coords[:, 0],
                UMAP_2=coords[:, 1],
            ),
        )


def add_multimodal_features(adata, mdata, features):
    if not adata.obs_names.isin(mdata.obs_names).all():
        raise ValueError("Clustered RNA cells absent from MuData")

    target = adata.obs_names

    for key in mdata.mod:
        flag = f"has_{key}"
        if flag in mdata.obs:
            features[flag] = (
                mdata.obs[flag].reindex(target).astype(bool)
                .astype(str).astype("category")
            )

    if "mod" in mdata.mod:
        obj = mdata.mod["mod"]
        depth = pd.Series(
            row_sums(obj.layers["callable_depth"]),
            index=obj.obs_names,
        )
        alt = pd.Series(
            row_sums(obj.layers["alt_count"]),
            index=obj.obs_names,
        )
        denominator = pd.Series(
            row_sums(obj.layers["ref_alt_depth"]),
            index=obj.obs_names,
        )

        features["mod_callable_depth"] = depth.reindex(target)
        features["mod_alt_fraction"] = (
            alt / denominator.where(denominator > 0)
        ).reindex(target)

    if "end5" in mdata.mod:
        obj = mdata.mod["end5"]
        features["end5_locus_counts"] = pd.Series(
            row_sums(obj.X),
            index=obj.obs_names,
        ).reindex(target)

        if "assigned_ends_fraction" in obj.obs:
            features["end5_assigned_fraction"] = (
                pd.to_numeric(
                    obj.obs["assigned_ends_fraction"],
                    errors="coerce",
                )
                .reindex(target)
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clustered", required=True)
    parser.add_argument("--mudata")
    parser.add_argument("--output", default="umap")
    parser.add_argument("--plot-format", default="png")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.clustered)

    if "X_umap" not in adata.obsm:
        raise ValueError("Clustered AnnData lacks X_umap")

    if not adata.obs_names.is_unique:
        raise ValueError("Duplicate clustered RNA barcodes")

    prefix = Path(args.output)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    coords = np.asarray(adata.obsm["X_umap"])
    groups = biotype_groups(adata)

    # Use sample-qualified observation IDs for multimodal joining.
    if "sample_id" not in adata.obs:
        raise ValueError("Missing RNA sample_id")

    if args.mudata:
        mdata = mu.read_h5mu(args.mudata)
        sample_ids = adata.obs["sample_id"].astype(str)
        qualified = pd.Index([
            f"{sample}:{barcode}"
            for sample, barcode in zip(sample_ids, adata.obs_names)
        ])
        adata.obs_names = qualified
    else:
        mdata = None

    obs_features = adata.obs.copy()

    general = [
        column for column in (
            "leiden", "total_counts", "n_genes_by_counts",
            "doublet_score", "predicted_doublet",
            "pct_counts_is_tRNA", "pct_counts_is_rRNA",
            "pct_counts_is_protein_coding", "tRNA_rRNA_ratio",
        ) if column in obs_features
    ]

    # Collect individual feature expression from log-normalised X.
    feature_expression = pd.DataFrame(index=adata.obs_names)

    for group in groups.values():
        for feature in group:
            j = adata.var_names.get_loc(feature)
            feature_expression[feature_label(adata, feature)] = (
                dense_column(adata.X, j)
            )

    if mdata is not None:
        add_multimodal_features(adata, mdata, obs_features)

    multi = [
        col for col in (
            "has_mod", "mod_callable_depth",
            "mod_alt_fraction", "has_end5",
            "end5_locus_counts", "end5_assigned_fraction",
        ) if col in obs_features
    ]

    plot_grid(
        coords, obs_features, general,
        f"{prefix}.general", args.plot_format,
    )

    for name, feature_ids in groups.items():
        labels = [
            feature_label(adata, feature)
            for feature in feature_ids
        ]
        plot_grid(
            coords, feature_expression, labels,
            f"{prefix}.{name}", args.plot_format,
        )

    plot_grid(
        coords, obs_features, multi,
        f"{prefix}.multimodal", args.plot_format,
    )

    print(
        f"Plotted {len(groups['trna_tmrna'])} tRNA/tmRNA "
        f"and {len(groups['ncrna'])} ncRNA features"
    )


if __name__ == "__main__":
    main()
