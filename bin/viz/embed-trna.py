
#!/usr/bin/env python3
"""Embed cells using relative tRNA expression alone."""

from typing import Callable
import argparse

import anndata as ad
from carabiner import print_err
from carabiner.mpl import grid, figsaver
import scanpy as sc
import numpy as np
import pandas as pd


def save_anndata(x, filename: str, printf: Callable = print) -> None:
    printf(f"Saving the following AnnData object as {filename}.h5ad:")
    printf(x)
    x.write(f"{filename}.h5ad", compression="gzip")
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--clustered", required=True)
    p.add_argument("--output", default="trna")
    p.add_argument("--min-cells-per-gene", type=int, default=3)
    p.add_argument("--plot-format", type=str, default="png")
    args = p.parse_args()

    rna = ad.read_h5ad(args.clustered)
    if "counts" not in rna.layers:
        raise ValueError("Missing raw counts")
    if "gene_biotype" not in rna.var:
        raise ValueError("Missing gene_biotype")
    if "is_tRNA" not in rna.var:
        raise ValueError("Missing is_tRNA")
    
    adata_trna = rna[:, rna.var["is_tRNA"].astype(bool)].copy()

    # Filter features using raw counts, not previously normalised X.
    adata_trna.X = adata_trna.layers["counts"].copy()
    sc.pp.filter_genes(
        adata_trna,
        min_cells=args.min_cells_per_gene,
    )

    if adata_trna.n_vars < 3:
        raise ValueError("Fewer than three qualifying tRNAs")

    depth = np.asarray(
        adata_trna.layers["counts"].sum(axis=1)
    ).ravel()

    observed = depth > 0
    if observed.sum() < 5:
        raise ValueError("Fewer than five tRNA-positive cells")

    n_excluded = int((~observed).sum())
    adata_trna = adata_trna[observed, :].copy()

    # Clear inherited log-normalisation metadata.
    adata_trna.uns.pop("log1p", None)

    sc.pp.normalize_total(adata_trna)
    sc.pp.log1p(adata_trna)

    print_err(f"tRNA PCA input: {adata_trna.shape}")

    n_pcs = min(30, adata_trna.n_obs - 1, adata_trna.n_vars - 1)
    if n_pcs < 2:
        raise ValueError("Insufficient tRNA PCA dimensions")

    from scipy.sparse import issparse

    if issparse(adata_trna.X) and adata_trna.n_vars < 100:
        adata_trna.X = adata_trna.X.toarray()
        assert adata_trna.X.shape == (adata_trna.n_obs, adata_trna.n_vars)
        assert np.isfinite(adata_trna.X).all()

    x = np.asarray(adata_trna.X)

    print_err(f"AnnData shape: {adata_trna.shape}")
    print_err(f"PCA matrix shape: {x.shape}")
    print_err(f"PCA matrix type: {type(adata_trna.X)}")
    print_err(f"PCA matrix dtype: {x.dtype}")
    print_err(f"tRNA matrix: {adata_trna.shape}")
    print_err(
        f"Inherited HVGs: "
        f"{adata_trna.var['highly_variable'].sum()}"
        if "highly_variable" in adata_trna.var
        else "No HVG annotation"
    )
    sc.pp.pca(
        adata_trna,
        n_comps=n_pcs,
        mask_var=None,
        svd_solver="full" if adata_trna.n_vars < 100 and not issparse(adata_trna.X) else "arpack",
        random_state=42,
    )
    sc.pp.neighbors(
        adata_trna,
        n_neighbors=min(15, adata_trna.n_obs - 1),
        use_rep="X_pca",
        random_state=42,
    )
    sc.tl.umap(adata_trna, random_state=42, min_dist=.1, spread=1.)

    adata_trna.obs["trna_depth"] = depth[observed]
    adata_trna.uns["trna_embedding"] = {
        "n_cells_input": rna.n_obs,
        "n_cells_excluded": n_excluded,
        "n_trna_features": adata_trna.n_vars,
        "normalization": "tRNA total, log1p",
    }

    save_anndata(
        adata_trna,
        f"{args.output}",
        printf=print_err,
    )

    data = adata_trna.obs.copy()
    data["UMAP_1"] = adata_trna.obsm["X_umap"][:, 0]
    data["UMAP_2"] = adata_trna.obsm["X_umap"][:, 1]

    colors = [
        c for c in (
            "leiden", "trna_depth", "total_counts",
            "pct_counts_is_tRNA", "n_genes_by_counts"
        ) if c in data
    ]

    fig, axes = grid(
        ncol=len(colors),
        aspect_ratio=1.,
        panel_size=3.75,
    )

    for ax, name in zip(fig.axes, colors):
        values = data[name]

        if isinstance(values.dtype, pd.CategoricalDtype):
            for category in values.cat.categories:
                selected = (values == category).to_numpy()
                ax.scatter(
                    *adata_trna.obsm["X_umap"][selected].T,
                    label=str(category),
                    s=1,
                    rasterized=True,
                )
            ax.legend(fontsize=7, markerscale=4)
        else:
            artist = ax.scatter(
                *adata_trna.obsm["X_umap"].T,
                c=values.to_numpy(),
                s=1,
                cmap="cividis",
                rasterized=True,
            )
            fig.colorbar(artist, ax=ax)

        ax.set(
            title=name,
            xlabel="tRNA UMAP 1",
            ylabel="tRNA UMAP 2",
        )

    figsaver(format=args.plot_format)(
        fig=fig,
        name=f"{args.output}.umap",
        df=pd.DataFrame(
            adata_trna.obsm['X_umap'], 
            columns=["UMAP_1", "UMAP_2"], 
            index=adata_trna.obs_names,
        ),
    )

    print_err(
        f"tRNA embedding: {adata_trna.n_obs} cells, "
        f"{adata_trna.n_vars} genes; excluded {n_excluded} cells"
    )


if __name__ == "__main__":
    main()
