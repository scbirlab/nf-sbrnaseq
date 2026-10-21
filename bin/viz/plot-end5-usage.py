
#!/usr/bin/env python3
"""Visualise alternative 5'-end locus usage."""

import argparse

import anndata as ad
from carabiner import print_err
from carabiner.mpl import grid, figsaver
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import seaborn as sns


def column(matrix, j):
    x = matrix[:, j]
    return x.toarray().ravel() if issparse(x) else np.asarray(x).ravel()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--clustered", required=True)
    p.add_argument("--mudata", required=True)
    p.add_argument("--output", default="end5")
    p.add_argument("--min-ends", type=int, default=3)
    p.add_argument("--top-loci", type=int, default=48)
    p.add_argument("--plot-format", type=str, default="png")
    args = p.parse_args()

    rna = ad.read_h5ad(args.clustered)
    mdata = mu.read_h5mu(args.mudata)

    # Match the existing human-readable MuData naming.
    mdata_end5 = mdata.mod["end5"]

    names = pd.Index([
        f"{sample}:{barcode}"
        for sample, barcode in zip(
            rna.obs["sample_id"].astype(str),
            rna.obs_names,
        )
    ])
    indices = mdata_end5.obs_names.get_indexer(names)
    if (indices < 0).any():
        raise ValueError("Missing RNA cells in 5'-end modality")

    counts = mdata_end5.X[indices]
    depth = np.asarray(counts.sum(axis=1)).ravel()
    present = mdata_end5.obs["present_in_modality"].astype(bool).values[indices]
    eligible = present & (depth >= args.min_ends)

    if not eligible.any() or mdata_end5.n_vars == 0:
        raise ValueError("No sufficiently covered 5'-end data")

    clusters = rna.obs["leiden"].astype(str).to_numpy()
    labels = sorted(np.unique(clusters))

    site_depth = np.asarray(counts[eligible].sum(axis=0)).ravel()
    order = np.argsort(-site_depth, kind="stable")
    selected = order[:args.top_loci]

    rows = []
    for j in selected:
        site = mdata_end5.var_names[j]
        these_counts = column(counts, j)
        record = {
            "locus_id": site,
            "total_count": int(these_counts.sum()),
        }

        for cluster in labels:
            mask = eligible & (clusters == cluster)
            denominator = depth[mask].sum()
            record[cluster] = (
                these_counts[mask].sum() / denominator
                if denominator else np.nan
            )

        rows.append(record)

    usage = pd.DataFrame(rows).set_index("locus_id")
    usage.to_csv(f"{args.output}.cluster-usage.csv")

    # Preserve the full locus annotation table.
    mdata_end5.var.to_csv(f"{args.output}.annotations.csv")

    def heatmap(frame, label, title, row_idx=None, col_idx=None):
        for_hm = frame.copy()
        if row_idx is not None:
            for_hm = for_hm.iloc[row_idx]
        if col_idx is not None:
            for_hm = for_hm.iloc[:, col_idx]

        missing = for_hm.isna()
        hm = sns.clustermap(
            for_hm.fillna(0.),
            metric="cosine",
            row_cluster=row_idx is None and len(for_hm) > 1,
            col_cluster=col_idx is None and len(for_hm.columns) > 1,
            cbar_kws={"label": label},
            cmap="cividis",
            vmin=0.,
            vmax=1.,
            xticklabels=True,
            yticklabels=True,
            mask=missing.to_numpy(),
            figsize=(
                max(6, min(14, for_hm.shape[1] * .4)),
                max(6, min(14, for_hm.shape[0] * .2)),
            ),
        )
        hm.ax_heatmap.set_yticklabels(
            hm.ax_heatmap.get_ymajorticklabels(), 
            fontsize=2.,
        )
        hm.savefig(
            f"{args.output}.{title}.{args.plot_format}", 
            dpi=300, 
            bbox_inches="tight",
        )
        return hm

    
    heatmap(
        usage,
        "Fraction of assigned ends",
        "Relative 5-end locus usage",
    )

    saver = figsaver(format=args.plot_format)
    xy = np.asarray(rna.obsm["X_umap"])

    # Site-specific UMAPs, paginated.
    for offset in range(0, len(selected), 24):
        page = selected[offset:offset + 24]
        fig, axes = grid(
            ncol=min(4, len(page)),
            nrow=int(np.ceil(len(page) / 4)),
            aspect_ratio=1.,
            panel_size=3.75,
        )

        data = pd.DataFrame(
            rna.obsm["X_umap"],
            columns=["UMAP_1", "UMAP_2"],
            index=names,
        )

        for ax, j in zip(fig.axes, page):
            site = mdata_end5.var_names[j]
            these_counts = column(counts, j)

            fraction = np.divide(
                these_counts,
                depth,
                out=np.full(len(these_counts), np.nan, dtype=float),
                where=eligible,
            )
            data[site] = fraction

            ax.scatter(
                *xy[~eligible].T,
                c="lightgrey",
                s=1,
                rasterized=True,
            )
            scatter = ax.scatter(
                *xy[eligible].T,
                c=fraction[eligible],
                cmap="cividis",
                vmin=0,
                vmax=1,
                s=1,
                rasterized=True,
            )
            fig.colorbar(scatter, ax=ax)
            ax.set(
                title=site,
                xlabel="UMAP 1",
                ylabel="UMAP 2",
            )

        saver(
            fig=fig,
            name=f"{args.output}.umap.page{offset // 24 + 1:02d}",
            df=data,
        )

    print_err(
        f"Plotted {len(selected)} loci; "
        f"{eligible.sum()} eligible cells"
    )


if __name__ == "__main__":
    main()
