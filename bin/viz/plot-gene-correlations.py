
#!/usr/bin/env python
"""Exploratory gene-gene correlation analysis."""

import argparse

import anndata as ad
from carabiner import print_err
from carabiner.mpl import grid, figsaver
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import seaborn as sns


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--clustered", required=True)
    p.add_argument("--output", default="correlations")
    p.add_argument("--min-cells-per-gene", type=int, default=3)
    p.add_argument("--max-genes", type=int, default=150)
    p.add_argument("--plot-format", type=str, default="png")
    args = p.parse_args()

    rna = ad.read_h5ad(args.clustered)

    counts = rna.layers["counts"]
    detected = (counts > 0).sum(axis=0).ravel()

    indices = np.flatnonzero(detected >= args.min_cells_per_gene)
    if len(indices) < 2:
        raise ValueError(f"Insufficient genes with at least {args.min_cells_per_gene=}")

    counts = rna.X[:, indices]
    if issparse(counts):
        counts = counts.toarray()
    else:
        counts = np.asarray(counts)

    variance = np.var(counts, axis=0)
    order = np.argsort(-variance, kind="stable")
    selected = order[:args.max_genes]

    frame = pd.DataFrame(
        counts[:, selected],
        index=rna.obs_names,
        columns=rna.var_names[indices[selected]],
    )

    # Constant genes have undefined correlations.
    frame = frame.loc[:, frame.nunique() > 1]
    if frame.shape[1] < 2:
        raise ValueError(f"Insufficient variable genes: {frame.shape=}")

    def heatmap(df, label):
        for_hm = df.fillna(0.)
        corr = for_hm.corr(method="spearman")
        distances = np.clip(
            1. - corr.fillna(0).to_numpy(), 0., 2.
        )
        np.fill_diagonal(distances, 0.)
        
        hm = sns.clustermap(
            corr,
            cbar_kws={"label": label},
            cmap="vlag",
            vmin=-1., 
            vmax=1.,
            xticklabels=True,
            yticklabels=True,
            dendrogram_ratio=(0., .1),
            figsize=(
                max(6, min(14, corr.shape[1] * .1)),
                max(6, min(14, corr.shape[0] * .1)), 
            ),
        )
        hm.ax_heatmap.set_yticklabels(
            hm.ax_heatmap.get_ymajorticklabels(), 
            fontsize=10.,
        )
        hm.savefig(
            f"correlations.{args.plot_format}", 
            dpi=300, 
            bbox_inches="tight",
        )
        return hm

    heatmap(
        frame,
        "Spearman $\\rho$",
    )
    print(f"Calculated correlations for {frame.shape[1]} genes")


if __name__ == "__main__":
    main()
