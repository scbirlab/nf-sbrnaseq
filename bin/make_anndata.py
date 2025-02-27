#!/usr/bin/env python

import sys

from functools import partial

import anndata as ad
from carabiner.mpl import add_legend, figsaver, grid, scattergrid
import pandas as pd
import numpy as np
import scanpy as sc

from scipy.sparse import csr_matrix
from pandas.api.types import CategoricalDtype

CELL_ID = "cell_barcode"
GENE_ID = "gene_id"
MIN_COUNTS_PER_CELL = 60
MIN_CELLS_PER_GENE = 3
MIN_GENES_TO_KEEP = 200
MIN_CELLS_TO_KEEP = 950

def make_anndata_from_table(counts_table: str) -> pd.DataFrame:
    print_err(f"Loading {counts_table} as a sparse matrix...")
    df = pd.read_csv(counts_table, sep="\t", low_memory=False)
    cell_cat = CategoricalDtype(sorted(df[CELL_ID].unique()), ordered=True)
    gene_cat = CategoricalDtype(sorted(df[GENE_ID].unique()), ordered=True)

    cell_idx = df[CELL_ID].astype(cell_cat).cat.codes
    gene_idx = df[GENE_ID].astype(gene_cat).cat.codes
    matrix = csr_matrix(
        (df["umi_count"].astype(float), (cell_idx, gene_idx)),
        shape=(cell_cat.categories.size, gene_cat.categories.size)
    )
    return df, pd.DataFrame.sparse.from_spmatrix(
        matrix,
        index=cell_cat.categories,
        columns=gene_cat.categories,
    )


def save_anndata(x, filename: str, printf: Callable = print) -> None:
    printf(f"Saving the following AnnData object as {filename}.h5ad:")
    printf(x)
    x.write(f"{filename}.h5ad", compression="gzip")
    return None


def main() -> None:

    sample_id, filename = sys.argv[:2]
    logfile = open(f"{sample_id}.scanpy.log", "w")
    print_err = partial(print, file=logfile)

    count_df, matrix_df = make_anndata_from_table(filename)
    adata = ad.AnnData(matrix_df)
    print_err(adata)
    adata.obs["sample_id"] = sample_id

    print_err("Annotating gene features...")
    gene_ann_columns = ["Chr", "locus_tag", "Name", "gene_biotype", "featurecounts_count", "Length"]
    gene_ann_df = count_df[[GENE_ID] + gene_ann_columns].drop_duplicates().set_index(GENE_ID)
    gene_ann_df = gene_ann_df.loc[adata.var_names,:]
    for col in gene_ann_columns:
        adata.var[col.casefold()] = gene_ann_df[col]

    all_biotypes = adata.var["gene_biotype"].unique()
    biotype_flags = []
    for biotype in all_biotypes:
        this_flag = f"is_{biotype}"
        biotype_flags.append(this_flag)
        adata.var[this_flag] = (adata.var["gene_biotype"] == biotype)
    print_err(adata)

    print_err("Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=biotype_flags, 
        inplace=True,
    )
    print_err(adata)

    print_err("Calculating tRNA:rRNA ratio...")
    adata.obs["tRNA_rRNA_ratio"] = (adata.obs["pct_counts_is_tRNA"] + 1) / (adata.obs["pct_counts_is_rRNA"] + 1)
    fig, axes = scattergrid(
        adata.obs,
        grid_columns=["n_genes_by_counts", "total_counts", "pct_counts_is_protein_coding", "pct_counts_is_tRNA", "pct_counts_is_rRNA", "tRNA_rRNA_ratio"],
        log=["n_genes_by_counts", "total_counts", ],
        aspect_ratio=1.25,
    )
    figsaver(format="png")(
        fig=fig,
        name='${sample_id}.cell-qc',
    )
    fig, axes = scattergrid(
        adata.var,
        grid_columns=["n_cells_by_counts", "total_counts", "length", "pct_dropout_by_counts",],
        log=["n_cells_by_counts", "total_counts", "length"],
        group="chr",
        aspect_ratio=1.25,
    )
    figsaver(format="png")(
        fig=fig,
        name='${sample_id}.gene-qc',
    )
    sc.pp.log1p(
        adata, 
    )

    save_anndata(adata, f"{sample_id}-unfiltered", printf=print_err)

    # Want to relax the cutoff with low counts to keep at least 1000 cells
    max_cutoff = adata.obs.nlargest(
        MIN_CELLS_TO_KEEP, 
        'total_counts',
    )['total_counts'].min()
    print_err(adata.obs.nlargest(
        MIN_CELLS_TO_KEEP, 
        'total_counts',
    ).sort_values('total_counts')['total_counts'])
    print_err(f"{max_cutoff=}")
    max_cutoff = max(1, min(int(max_cutoff), MIN_COUNTS_PER_CELL))
    print_err(f"{max_cutoff=}")
    print_err(f"Filtering cells with {max_cutoff=}...")
    adata = adata[adata.obs['total_counts'] >= max_cutoff]

    print_err(adata)
    # Want to relax the cutoff with low counts to keep at least {MIN_GENES_TO_KEEP} genes
    max_cutoff = adata.var.nlargest(
        MIN_GENES_TO_KEEP + 1, 
        'n_cells_by_counts',
    )['n_cells_by_counts'].min()
    print_err(f"{max_cutoff=}")
    max_cutoff = max(1, min(int(max_cutoff) - 1, MIN_CELLS_PER_GENE))
    print_err(f"{max_cutoff=}")
    print_err(f"Filtering genes with {max_cutoff=}...")
    sc.pp.filter_genes(
        adata, 
        min_cells=max_cutoff,
    )
    print_err(adata)

    #print_err(f"Removing rRNA...")
    #adata = adata[:, ~adata.var["is_rRNA"]]
    #print_err(adata)
    print_err("Normalizing to counts per cell...")
    sc.pp.normalize_total(
        adata, 
        exclude_highly_expressed=True,
    )
    print_err(f"Annotating highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat',
    )
    adata.var["not_is_rRNA"] = (~adata.var["is_rRNA"])
    adata.var["highly_variable_non_rRNA"] = ((~adata.var["is_rRNA"]) & adata.var["highly_variable"])
    print_err(f"Found {adata.var['highly_variable_non_rRNA'].sum()=} highly variable genes!")
    print_err(adata)
    use_highly_variable = (adata.var['highly_variable_non_rRNA'].sum() >= 3)

    print_err("Applying PCA...")
    sc.pp.pca(
        adata, 
        mask_var="highly_variable_non_rRNA" if use_highly_variable else "not_is_rRNA",
    )
    print_err("Getting nearest neighbours...")
    sc.pp.neighbors(
        adata, 
        metric='cosine',
    )
    print_err("Carrying out Leiden clustering on cells...")
    sc.tl.leiden(
        adata, 
        key_added="leiden", 
        resolution=1.,
        flavor='igraph',
        directed=False,
    )
    print_err(f"Found {adata.obs['leiden'].cat.categories.size=} Leiden clusters of cells!")
    print_err("Finding gene markers per cell cluster...")
    sc.tl.rank_genes_groups(
        adata, 
        "leiden",
        mask_var="highly_variable_non_rRNA" if use_highly_variable else "not_is_rRNA",
        method="logreg",
    )
    print_err("Embedding as UMAP...")
    sc.tl.umap(
        adata, 
        min_dist=.5,  # umap-learn default = .1
        spread=1.,
    )

    save_anndata(adata, sample_id, printf=print_err)

    n_genes = adata.n_vars
    n_cells = adata.n_obs

    colors_to_plot = [
        "total_counts", 
        "n_genes_by_counts",
        "doublet_score",
        "predicted_doublet",
        "leiden",
        "tRNA_rRNA_ratio",
    ] + [f"pct_counts_{biotype}" for biotype in biotype_flags]
    colors_to_plot = [c for c in colors_to_plot if c in adata.obs]
    log_colors = [c for c in colors_to_plot if c.endswith("_counts")] + ["tRNA_rRNA_ratio",]
    print_err("Making UMAP plots...")
    fig, axes = grid(ncol=len(colors_to_plot), aspect_ratio=1., panel_size=3.75)
    for ax, col in zip(axes, colors_to_plot):
        colors = adata.obs[col]
        col_is_categorical = (isinstance(colors, pd.Series) and colors.dtype == "category")
        use_log = (not col_is_categorical and np.any(colors > 0.) and col in log_colors)
        plotter = partial(
            ax.scatter,
            s=.5,
            alpha=.7,
            plotnonfinite=True,
        )
        if col_is_categorical:
            for i in np.unique(colors.cat.codes):
            plotter(
                *adata[adata.obs[col].cat.codes == i].obsm['X_umap'].T,
                label=i,
            )
            ax.legend(
            loc='upper center',
            bbox_to_anchor=(0.5, -0.2),
            ncol=7,

            )
        if not col_is_categorical:
            scatter = plotter(
            *adata.obsm['X_umap'].T,
            c=colors,
            cmap="cividis",
            norm="log" if use_log else None,
            )
            fig.colorbar(scatter, ax=ax)
        ax.set(
            title=f"{n_cells} cells x {n_genes} genes\\n{col}", 
            xlabel="UMAP_1", 
            ylabel="UMAP_2",
        )
    print_err(f"Saving UMAP plots as {sample_id}.umap...")
    figsaver(format="png")(
        fig=fig,
        name=f'{sample_id}.umap',
    )
    print_err("Plotting heatmap of mean gene markers per cell cluster...")
    sc.pl.rank_genes_groups_matrixplot(
        adata,
        save=f"{sample_id}.cluster-gene-markers-mean.png",
    )
    print_err("Plotting heatmap of gene markers per cell cluster...")
    sc.pl.rank_genes_groups_heatmap(
        adata, 
        show_gene_labels=True, 
        save=f"{sample_id}.cluster-gene-markers.png",
    )
    print_err("Done!")
    return None

if __name__ == '__main__':
    main()