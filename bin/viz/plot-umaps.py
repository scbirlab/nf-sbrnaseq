
#!/usr/bin/env python
"""Plot RNA and multimodal UMAP feature grids."""

import argparse
from pathlib import Path

import anndata as ad
from carabiner.mpl import add_legend, grid, figsaver
from functools import partial
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

    biotype = adata.var["gene_biotype"].str.casefold()

    trna = biotype.eq("trna")
    tmrna = biotype.isin(["tmrna", "transfer-messenger rna"])
    rrna = biotype.eq("rrna")

    other_ncrna = (
        biotype.str.contains(
            "|".join([
                "ncrna", "non_coding", "non-coding", "srna", "antisense"
        ]))
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
            aspect_ratio=1.,
            panel_size=3.75,
        )

        for ax, col in zip(fig.axes, page):
            colors = frame[col]
            valid = colors.notna().to_numpy()
            col_is_categorical = (
                isinstance(colors, pd.Series) 
                and colors.dtype in ("category", "object", "str")
            )
            use_log = (not col_is_categorical and np.all(colors > 0.))
            plotter = partial(
                ax.scatter,
                s=.5,
                plotnonfinite=True,
            )
            if col_is_categorical:
                colors = colors.astype("category")
                unique_cats = np.unique(colors)
                min_cat = sorted(colors)[0]
                for i in np.unique(colors):
                    if len(unique_cats) == 2 and i == min_cat:
                        plotter(
                            *coords[colors.cat.codes == i].T,
                            label=i,
                            color="lightgrey",
                        )
                    plotter(
                        *coords[colors.cat.codes == i].T,
                        label=i,
                    )
                add_legend(ax)
            if not col_is_categorical:
                mask = np.isfinite(colors)
                if np.all(colors >= 0.):
                    mask = (colors != 0.) & mask
                scatter = plotter(
                    *coords[mask].T,
                    c=colors[mask],
                    cmap="cividis",
                    norm="log" if use_log else None,
                )
                fig.colorbar(scatter, ax=ax)

            ax.set(
                title=f"{coords.shape[0]} cells" + '\n' + col,
                xlabel="UMAP 1",
                ylabel="UMAP 2",
            )

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
    _biotype_groups = biotype_groups(adata)

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

    general = [
        column for column in (
            "leiden", "total_counts", "n_genes_by_counts",
            "doublet_score", "predicted_doublet",
            "pct_counts_is_tRNA", "pct_counts_is_rRNA",
            "pct_counts_is_protein_coding", "tRNA_rRNA_ratio",
        ) if column in adata.obs
    ]

    # Collect individual feature expression from log-normalised X.
    feature_expression = pd.DataFrame(index=adata.obs_names)

    for feature_ids in _biotype_groups.values():
        for feature in feature_ids:
            j = adata.var_names.get_loc(feature)
            feature_expression[feature_label(adata, feature)] = (
                dense_column(adata.X, j)
            )

    if mdata is not None:
        add_multimodal_features(adata, mdata, adata.obs)

    multi = [
        col for col in (
            "has_mod", "mod_callable_depth",
            "mod_alt_fraction", "has_end5",
            "end5_locus_counts", "end5_assigned_fraction",
        ) if col in adata.obs
    ]

    plot_grid(
        coords, 
        adata.obs, 
        general,
        f"{prefix}.general", 
        args.plot_format,
    )

    for name, feature_ids in _biotype_groups.items():
        labels = [
            feature_label(adata, feature)
            for feature in feature_ids
        ]
        plot_grid(
            coords, 
            feature_expression, 
            labels,
            f"{prefix}.{name}", 
            args.plot_format,
        )

    plot_grid(
        coords, 
        adata.obs, 
        multi,
        f"{prefix}.multimodal", 
        args.plot_format,
    )

    print(
        f"Plotted {len(_biotype_groups['trna_tmrna'])} tRNA/tmRNA "
        f"and {len(_biotype_groups['ncrna'])} ncRNA features"
    )


if __name__ == "__main__":
    main()
