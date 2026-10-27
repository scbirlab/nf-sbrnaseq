
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

from viz_utils import heatmap, end_locus_label


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

    # Always preserve locus annotations, even when plotting is skipped.
    mdata_end5.var.to_csv(f"{args.output}.annotations.csv")

    reason = None
    if mdata_end5.n_vars == 0:
        reason = "No 5'-end loci detected"
    elif not eligible.any():
        reason = (
            f"No cells with at least {args.min_ends} "
            "assigned 5'-end molecules"
        )

    pd.DataFrame([{
        "status": "skipped" if reason else "completed",
        "reason": reason or "",
        "n_rna_cells": rna.n_obs,
        "n_loci": mdata_end5.n_vars,
        "n_cells_present": int(present.sum()),
        "n_cells_eligible": int(eligible.sum()),
        "total_assigned_ends": int(depth.sum()),
        "min_ends": args.min_ends,
    }]).to_csv(f"{args.output}.status.csv", index=False)

    if reason:
        print_err(f"Skipping 5'-end visualisation: {reason}")
        return None

    clusters = rna.obs["leiden"].astype(str).to_numpy()
    labels = sorted(np.unique(clusters))

    
    site_depth = np.asarray(counts[eligible].sum(axis=0)).ravel()

    # Select the most abundant loci first.
    selected = np.argsort(
        -site_depth, kind="stable"
    )[:args.top_loci]

    # Display selected loci in genomic coordinate order.
    positions = mdata_end5.var.iloc[selected].copy()
    positions["start1"] = pd.to_numeric(
        positions["start1"], errors="raise"
    )
    positions["end1"] = pd.to_numeric(
        positions["end1"], errors="raise"
    )
    positions["_index"] = selected

    positions = positions.sort_values(
        ["contig", "start1", "end1"],
        kind="stable",
    )
    selected = positions["_index"].to_numpy(dtype=int)

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

    # Display names are separate from the stable locus IDs.
    display_labels = [
        end_locus_label(
            site,
            mdata_end5.var.loc[site].get(
                "feature_attributes", ""
            ),
        )
        for site in usage.index
    ]

    display_usage = usage[labels].copy()
    display_usage.index = display_labels

    # Rows are already in genomic coordinate order.
    # Disable row clustering; retain column clustering.
    hm = heatmap(
        display_usage,
        color_label="Fraction of assigned ends",
        filename=f"{args.output}.end5-locus-usage",
        output_format=args.plot_format,
        row_idx=np.arange(len(display_usage)),
        figsize=(
            max(8, len(labels) * 1.1),
            max(8, len(display_usage) * .35),
        ),
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
                title=end_locus_label(
                    site,
                    mdata_end5.var.loc[site].get("feature_attributes", ""),
                ),
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
