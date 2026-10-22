
#!/usr/bin/env python
"""Modification-site UMAP grids and cluster-level heatmaps."""

import argparse
from pathlib import Path

import anndata as ad
from carabiner import print_err
from carabiner.mpl import grid, figsaver
from matplotlib.colors import Normalize
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import seaborn as sns


def vector(matrix, rows, column):
    """Read one feature without densifying the full matrix."""
    x = matrix[rows, column]
    return x.toarray().ravel() if issparse(x) else np.asarray(x).ravel()


def run(clustered, mudata, prefix, fmt, min_depth, page_size):
    rna = ad.read_h5ad(clustered)
    mdata = mu.read_h5mu(mudata)

    if "X_umap" not in rna.obsm or "leiden" not in rna.obs:
        raise ValueError("RNA clustering or UMAP missing")
    if "mod" not in mdata.mod:
        raise ValueError("Modification modality missing")

    mod = mdata.mod["mod"]
    for name in (
        "ref_count", 
        "alt_count",
        "callable_depth", 
        "ref_alt_depth",
    ):
        if name not in mod.layers:
            raise ValueError(f"Missing modification layer: {name}")

    if "sample_id" not in rna.obs:
        raise ValueError("RNA sample_id missing")
    if not rna.obs_names.is_unique or not mod.obs_names.is_unique:
        raise ValueError("Duplicate barcode identifiers")

    names = pd.Index([
        f"{sample}:{barcode}"
        for sample, barcode in zip(
            rna.obs["sample_id"].astype(str),
            rna.obs_names,
        )
    ])

    indices = mod.obs_names.get_indexer(names)
    if np.any(indices < 0):
        raise ValueError("RNA cells absent from MuData")

    present = mod.obs["present_in_modality"].to_numpy(dtype=bool)[indices]

    coords = np.asarray(rna.obsm["X_umap"])
    clusters = rna.obs["leiden"].astype(str).to_numpy()
    cluster_names = sorted(
        np.unique(clusters),
        key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x),
    )

    saver = figsaver(format=fmt)
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    sites = mod.var_names.tolist()
    summary = []
    coverage_rows = []
    fraction_rows = []

    for j, site in enumerate(sites):
        ref = vector(mod.layers["ref_count"], indices, j)
        alt = vector(mod.layers["alt_count"], indices, j)
        depth = vector(mod.layers["callable_depth"], indices, j)
        denominator = ref + alt

        frac = np.divide(
            alt, denominator,
            out=np.full(len(alt), np.nan, dtype=float),
            where=(denominator >= min_depth) & present,
        )

        row = {
            "site_id": site,
            "n_cells_callable": int(np.sum((depth > 0) & present)),
            "n_cells_ref_alt": int(
                np.sum((denominator >= min_depth) & present)
            ),
            "total_callable_depth": int(depth.sum()),
            "total_ref_count": int(ref.sum()),
            "total_alt_count": int(alt.sum()),
        }
        for field in (
            "feature_id", "feature", "ref_base",
            "alt_base", "modification",
        ):
            if field in mod.var:
                row[field] = mod.var.iloc[j][field]
        summary.append(row)

        coverage_row = {"site_id": site}
        fraction_row = {"site_id": site}

        for cluster in cluster_names:
            selected = (clusters == cluster) & present
            n = int(selected.sum())
            coverage_row[cluster] = (
                float(np.mean(depth[selected] > 0)) if n else np.nan
            )
            denom_sum = denominator[selected].sum()
            fraction_row[cluster] = (
                float(alt[selected].sum() / denom_sum)
                if denom_sum >= min_depth else np.nan
            )

        coverage_rows.append(coverage_row)
        fraction_rows.append(fraction_row)

    if not sites:
        print_err("No modification sites detected; analysis completed.")
        return None

    pd.DataFrame(summary).to_csv(
        f"{prefix}.sites.csv", index=False
    )

    coverage_df = pd.DataFrame(coverage_rows).set_index("site_id")
    fraction_df = pd.DataFrame(fraction_rows).set_index("site_id")

    coverage_df.to_csv(f"{prefix}.cluster-coverage.csv")
    fraction_df.to_csv(f"{prefix}.cluster-alt-fraction.csv")


    def heatmap(df, label, title, row_idx=None, col_idx=None):
        if row_idx is not None:
            df = df.iloc[row_idx]
        if col_idx is not None:
            df = df.iloc[:,col_idx]
        missing = df.isna()
        values = df.fillna(0.)
        # Only cluster an axis containing distinguishable observations.
        row_cluster = (
            row_idx is None
            and values.shape[0] > 1
            and values.drop_duplicates().shape[0] > 1
        )
        col_cluster = (
            col_idx is None
            and values.shape[1] > 1
            and values.T.drop_duplicates().shape[0] > 1
        )
        hm = sns.clustermap(
            values,
            metric="cosine",
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            cbar_kws={"label": label},
            cmap="cividis",
            vmin=0., 
            vmax=1.,
            xticklabels=True,
            yticklabels=True,
            mask=missing,
            figsize=(
                max(6, min(14, values.shape[1] * .5)),
                max(6, min(14, values.shape[0] * .25)),
            ),
        )
        hm.ax_heatmap.tick_params(axis="y", labelsize=6)
        hm.savefig(
            f"{prefix}.{title}.png", 
            dpi=300, 
            bbox_inches="tight",
        )
        hm.savefig(
            f"{prefix}.{title}.pdf", 
            dpi=300, 
            bbox_inches="tight",
        )
        return hm

    hm = heatmap(
        fraction_df,
        "Pooled ALT / (REF + ALT)",
        "Site ALT fraction by cluster",
    )
    row_order = (
        hm.dendrogram_row.reordered_ind
        if hm.dendrogram_row is not None
        else None
    )
    col_order = (
        hm.dendrogram_col.reordered_ind
        if hm.dendrogram_col is not None
        else None
    )
    heatmap(
        coverage_df,
        "Fraction callable",
        "Site coverage by cluster",
        row_idx=hm.dendrogram_row.reordered_ind,
        col_idx=hm.dendrogram_col.reordered_ind,
    )

    # Each feature gets one coverage panel and one ALT-fraction panel.
    def umap_pages(kind):
        for offset in range(0, len(sites), page_size):
            page = sites[offset:offset + page_size]
            ncol = min(6, len(page))
            fig, axes = grid(
                ncol=ncol,
                nrow=int(np.ceil(len(page) / ncol)),
                aspect_ratio=1.1,
                panel_size=3.75,
            )

            plot_data = pd.DataFrame({
                "cell_id": names,
                "UMAP_1": coords[:, 0],
                "UMAP_2": coords[:, 1],
                "leiden": clusters,
            })

            for ax, site in zip(fig.axes, page):
                j = mod.var_names.get_loc(site)

                depth = vector(
                    mod.layers["callable_depth"], indices, j
                )
                ref = vector(mod.layers["ref_count"], indices, j)
                alt = vector(mod.layers["alt_count"], indices, j)
                denominator = ref + alt

                if kind == "coverage":
                    values = np.where(present, depth, np.nan)
                    valid = present & (depth > 0)
                    norm = "log"
                else:
                    valid = present & (denominator >= min_depth)
                    values = np.divide(
                        alt, denominator,
                        out=np.full(len(alt), np.nan),
                        where=valid,
                    )
                    if kind == "fraction":
                        norm = Normalize(
                            vmin=0,
                            vmax=1 if kind == "fraction"
                            else max(1, float(np.nanmax(values[valid]))),
                        )
                color = values

                plot_data[site] = values

                # Unmeasured or insufficiently covered cells.
                ax.scatter(
                    *coords[~valid].T,
                    s=.5, 
                    c="lightgrey",
                    rasterized=True,
                    zorder=-1,
                )

                if np.any(valid):
                    sc = ax.scatter(
                        *coords[valid].T,
                        c=color[valid],
                        cmap="cividis",
                        norm=norm,
                        s=.5,
                        rasterized=True,
                        zorder=1,
                    )
                    fig.colorbar(sc, ax=ax)

                ax.set(
                    title=str(site),
                    xlabel="UMAP 1",
                    ylabel="UMAP 2",
                )

            saver(
                fig=fig,
                name=(
                    f"{prefix}.umap-{kind}."
                    f"page{offset // page_size + 1:02d}"
                ),
                df=plot_data,
            )

    umap_pages("coverage")
    umap_pages("fraction")

    print_err(
        f"Plotted {len(sites)} modification sites "
        f"across {len(names)} RNA-QC-passing cells"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clustered", required=True)
    parser.add_argument("--mudata", required=True)
    parser.add_argument("--output", default="modsites")
    parser.add_argument("--plot-format", default="png")
    parser.add_argument("--min-depth", type=int, default=1)
    parser.add_argument("--page-size", type=int, default=24)
    args = parser.parse_args()

    if args.min_depth < 1 or args.page_size < 1:
        parser.error("min-depth and page-size must be positive")

    run(
        args.clustered, args.mudata, args.output,
        args.plot_format, args.min_depth, args.page_size,
    )


if __name__ == "__main__":
    main()
