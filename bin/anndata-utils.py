#!/usr/bin/env python
from typing import Callable

from argparse import Namespace, FileType
from functools import partial
import sys

import anndata as ad
from carabiner import print_err
from carabiner.cliutils import CLIApp, CLICommand, CLIOption, clicommand
from carabiner.mpl import add_legend, figsaver, grid, scattergrid
import pandas as pd
import numpy as np
import scanpy as sc

from scipy.sparse import csr_matrix, coo_matrix
from pandas.api.types import CategoricalDtype

CELL_ID = ["cell_barcode"]
GENE_ID = ["genome_accession", "gene_id"]
UMI_COUNT_COLUMN = "umi_count"
MIN_COUNTS_PER_CELL = 60
MIN_CELLS_PER_GENE = 3
MIN_GENES_TO_KEEP = 200
MIN_CELLS_TO_KEEP = 950

__version__ = "0.0.1"

def make_anndata_from_table(counts_table: str) -> pd.DataFrame:

    def _joiner(col: str, char: str = "-"):
        def f(df: pd.DataFrame) -> pd.Series:
            return df[col].apply(
                char.join, 
                raw=True, 
                axis=1, 
                result_type="reduce",
            )
        return f
    
    df = (
        pd.read_csv(
            counts_table, 
            sep="\t", 
            low_memory=False,
        )
        .assign(**{
            "__cell_id__": _joiner(CELL_ID),
            "__gene_id__": _joiner(GENE_ID),
        })
    )
    
    print_err(df)
    
    cell_cat = CategoricalDtype(sorted(df["__cell_id__"].unique()), ordered=True)
    gene_cat = CategoricalDtype(sorted(df["__gene_id__"].unique()), ordered=True)

    cell_idx = df["__cell_id__"].astype(cell_cat).cat.codes
    gene_idx = df["__gene_id__"].astype(gene_cat).cat.codes
    matrix = coo_matrix(
        (
            df[UMI_COUNT_COLUMN].astype(float), 
            (cell_idx, gene_idx),
        ),
        shape=(cell_cat.categories.size, gene_cat.categories.size),
    ).tocsr()
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


def plot_barnyard(
    adata,
    filename: str,
    grouping="genome_accession",
    plot_format: str = "png"
):
    _figsaver = figsaver(format=plot_format)
    # If more than 1 assembly, group by and plot scatter
    unique_assembly = sorted(adata.var[grouping].unique())
    n_assembly = len(unique_assembly)
    print_err(f"Found {n_assembly} unique genome reference(s): {unique_assembly}")
    if n_assembly > 1:
        assembly_counts = sc.get.aggregate(
            adata, 
            by=grouping, 
            func="sum", 
            axis="var",
        )
        assembly_counts = assembly_counts.to_df(layer="sum")
        print_err(assembly_counts)
        fig, axes = scattergrid(
            assembly_counts,
            grid_columns=assembly_counts.columns,
            aspect_ratio=1.25,
        )
        print_err("Saving barnyard plots...")
        _figsaver(
            fig=fig,
            name=filename, 
            df=assembly_counts.reset_index(),
        )
    return None


@clicommand(message="Building AnnData with the following parameters")
def _build(args: Namespace) -> None:

    _figsaver = figsaver(format=args.plot_format)
    print_err(f"Loading file {args.input} as a sparse matrix...")
    count_df, matrix_df = make_anndata_from_table(args.input)

    adata = ad.AnnData(matrix_df)
    adata.obs_names = matrix_df.index
    adata.var_names = matrix_df.columns
    adata.obs["sample_id"] = args.id

    print_err("Annotating gene features...")
    gene_ann_columns = GENE_ID + [
        "Chr", 
        "locus_tag", 
        "Name", 
        "gene_biotype", 
        "bulk_read_count", 
        "Length",
    ]
    gene_ann_df = (
        count_df[["__gene_id__"] + gene_ann_columns]
        .drop_duplicates()
        .set_index("__gene_id__")
    )
    gene_ann_df = gene_ann_df.loc[adata.var_names,:]
    for col in gene_ann_columns:
        adata.var[col.casefold()] = gene_ann_df[col]

    all_biotypes = adata.var["gene_biotype"].unique()
    biotype_flags = []
    for biotype in all_biotypes:
        this_flag = f"is_{biotype}"
        biotype_flags.append(this_flag)
        adata.var[this_flag] = adata.var["gene_biotype"] == biotype

    plot_barnyard(
        adata,
        filename=f"{args.output}.genome-counts",
        grouping="genome_accession",
        plot_format=args.plot_format,
    )

    print_err("Calculating QC metrics...")
    percentiles = [q for q in (50, 100) if q <= adata.n_vars]
    if not percentiles and adata.n_vars > 0:
        percentiles = [min(10, adata.n_vars)]
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=biotype_flags, 
        percent_top=percentiles,
        inplace=True,
    )

    print_err("Calculating tRNA:rRNA ratio...")
    adata.obs["tRNA_rRNA_ratio"] = (
        (adata.obs["pct_counts_is_tRNA"] + 1) 
        / 
        (adata.obs["pct_counts_is_rRNA"] + 1)
    )
    fig, axes = scattergrid(
        adata.obs,
        grid_columns=[
            "n_genes_by_counts", 
            "total_counts", 
            "pct_counts_is_protein_coding", 
            "pct_counts_is_tRNA", 
            "pct_counts_is_rRNA", 
            "tRNA_rRNA_ratio",
        ],
        log=[
            "n_genes_by_counts", 
            "total_counts",
        ],
        aspect_ratio=1.25,
    )
    _figsaver(
        fig=fig,
        name=f"{args.output}.cell-qc",
        df=adata.obs,
    )

    fig, axes = scattergrid(
        adata.var,
        grid_columns=[
            "n_cells_by_counts", 
            "total_counts", 
            "length", 
            "pct_dropout_by_counts",
        ],
        log=[
            "n_cells_by_counts", 
            "total_counts",
        ],
        grouping="genome_accession",
        aspect_ratio=1.25,
    )
    _figsaver(
        fig=fig,
        name=f"{args.output}.gene-qc",
        df=adata.var,
    )
    sc.pp.log1p(adata)

    save_anndata(adata, f"{args.output}", printf=print_err)

    return None

    
@clicommand(message="Filtering AnnData with the following parameters")
def _filter(args: Namespace) -> None:

    adata = sc.read_h5ad(args.input)

    # Want to relax the cutoff with low counts to keep at least 1000 cells
    max_cutoff = adata.obs.nlargest(
        args.min_cells, 
        'total_counts',
    )['total_counts'].min()
    max_cutoff = max(1, min(int(max_cutoff), args.min_counts_per_cell))
    print_err(f"Filtering cells with {max_cutoff=}...")
    adata = adata[adata.obs['total_counts'] >= max_cutoff]

    # Want to relax the cutoff with low counts to keep at least {MIN_GENES_TO_KEEP} genes
    max_cutoff = adata.var.nlargest(
        args.min_genes + 1, 
        'n_cells_by_counts',
    )['n_cells_by_counts'].min()
    max_cutoff = max(1, min(int(max_cutoff) - 1, args.min_cells_per_gene))
    print_err(f"Filtering genes with {max_cutoff=}...")
    sc.pp.filter_genes(
        adata, 
        min_cells=max_cutoff,
    )

    save_anndata(adata, f"{args.output}.filtered", printf=print_err)

    plot_barnyard(
        adata,
        filename=f"{args.output}.genome-counts-after-filter",
        grouping="genome_accession",
        plot_format=args.plot_format,
    )

    return None


@clicommand(message="Filtering AnnData with the following parameters")
def _cluster(args: Namespace) -> None:

    _figsaver = figsaver(format=args.plot_format)
    adata = sc.read_h5ad(args.input)
    print_err(adata)

    print_err("Annotating highly variable genes...")
    adata.var["not_is_rRNA"] = (~adata.var["is_rRNA"])
    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor='seurat',
        )
    except ValueError as e:
        print_err(e)
        use_highly_variable = False
    else:
        adata.var["highly_variable_non_rRNA"] = (adata.var["not_is_rRNA"] & adata.var["highly_variable"])
        print_err(f"Found {adata.var['highly_variable_non_rRNA'].sum()} highly variable genes!")
        use_highly_variable = adata.var['highly_variable_non_rRNA'].sum() >= 3

    print_err("Applying PCA...")
    try:
        sc.pp.pca(
            adata, 
            mask_var="highly_variable_non_rRNA" if use_highly_variable else "not_is_rRNA",
        )
    except ValueError as e:
        print_err(e)
    print_err("Getting nearest neighbours...")
    sc.pp.neighbors(
        adata, 
        metric='euclidean',
    )
    print_err("Carrying out Leiden clustering on cells...")
    sc.tl.leiden(
        adata, 
        key_added="leiden", 
        resolution=1.,
        flavor='igraph',
        directed=False,
    )
    print_err(f"Found {adata.obs['leiden'].cat.categories.size} Leiden clusters of cells!")
    print_err("Finding gene markers per cell cluster...")
    sc.tl.rank_genes_groups(
        adata, 
        "leiden",
        mask_var="highly_variable_non_rRNA" if use_highly_variable else "not_is_rRNA",
        method="t-test",
        rankby_abs=True,
    )
    print_err("Embedding as UMAP...")
    sc.tl.umap(
        adata, 
        min_dist=.1,  # umap-learn default = .1
        spread=1.,
    )
    save_anndata(adata, f"{args.output}.clustered", printf=print_err)

    n_genes, n_cells = adata.n_vars, adata.n_obs
    biotype_flags = [c for c in adata.obs if c.startswith("is_")]
    colors_to_plot = [
        "total_counts", 
        "n_genes_by_counts",
        "doublet_score",
        "predicted_doublet",
        "leiden",
        "tRNA_rRNA_ratio",
    ] + [
        c for c in adata.obs if c.startswith("pct_counts_")
    ]
    colors_to_plot = [c for c in colors_to_plot if c in adata.obs]
    log_colors = [c for c in colors_to_plot if c.endswith("_counts")] + ["tRNA_rRNA_ratio"]

    print_err("Making UMAP plots...")
    fig, axes = grid(
        ncol=len(colors_to_plot), 
        aspect_ratio=1., 
        panel_size=3.75,
    
    )
    for ax, col in zip(axes, colors_to_plot):
        colors = adata.obs[col]
        col_is_categorical = isinstance(colors, pd.Series) and colors.dtype in ("category", "object", "str")
        use_log = (not col_is_categorical and np.any(colors > 0.) and col in log_colors)
        plotter = partial(
            ax.scatter,
            s=.5,
            alpha=.7,
            plotnonfinite=True,
        )
        if col_is_categorical:
            colors = colors.astype("category")
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
            title=f"{n_cells} cells x {n_genes} genes" + '\n' + col, 
            xlabel="UMAP_1", 
            ylabel="UMAP_2",
        )
    _figsaver(
        fig=fig,
        name=f'{args.output}.umap',
        df=pd.DataFrame(adata.obsm['X_umap'], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names),
    )
    print_err("Plotting heatmap of mean gene markers per cell cluster...")
    # Compute a dendrogram on a dense representation only (avoids pandas+csr bug)
    try:
        if "X_pca" in adata.obsm_keys():
            sc.tl.dendrogram(adata, groupby="leiden", use_rep="X_pca", key_added="dendrogram_leiden")
            dendro_arg = "dendrogram_leiden"
        else:
            if sp.issparse(adata.X):
                adata.obsm["X_dense"] = adata.X.toarray()
                rep = "X_dense"
            else:
                rep = "X"
            sc.tl.dendrogram(adata, groupby="leiden", use_rep=rep, key_added="dendrogram_leiden")
            dendro_arg = "dendrogram_leiden"
    except Exception as e:
        print(f"[warn] dendrogram computation failed: {e}. Proceeding without dendrogram.")
        dendro_arg = False
    sc.pl.rank_genes_groups_matrixplot(
        adata,
        n_genes=min(10, adata.n_vars),
        groupby="leiden",
        dendrogram=dendro_arg,
        save=f"{args.output}.cluster-gene-markers-mean.{args.plot_format}",
    )
    print_err("Plotting heatmap of gene markers per cell cluster...")
    sc.pl.rank_genes_groups_heatmap(
        adata, 
        show_gene_labels=True, 
        save=f"{args.output}.cluster-gene-markers.{args.plot_format}",
    )

    return None


def main() -> None:

    min_cells = CLIOption(
        '--min-cells', '-c', 
        type=int, 
        default=MIN_CELLS_TO_KEEP,
        help='Minimum number of cells to retain after filtering.',
    )
    min_genes = CLIOption(
        '--min-genes', '-g', 
        type=int, 
        default=MIN_GENES_TO_KEEP,
        help='Minimum number of genes to retain after filtering.',
    )
    min_counts_per_cell = CLIOption(
        '--min-counts-per-cell', '-C', 
        type=int, 
        default=MIN_COUNTS_PER_CELL,
        help='Minimum number of counts per cell.',
    )
    min_cells_per_gene = CLIOption(
        '--min-cells-per-gene', '-G', 
        type=int, 
        default=MIN_CELLS_PER_GENE,
        help='Minimum number of cells per gene.',
    )

    inputs = CLIOption(
        'input', 
        type=str,
        nargs='?',
        help='Input file.',
    )
    sample_id = CLIOption(
        '--id', '-x', 
        type=str,
        help='Sample ID to annotate.',
    )
    plot_format = CLIOption(
        '--plot-format', '-p', 
        type=str,
        choices=["png", "pdf"],
        default="png",
        help='File type for plots.',
    )

    outputs = CLIOption(
        '--output', '-o', 
        type=str,
        default="output",
        help='Output file prefix.',
    )

    build = CLICommand(
        "build", 
        description="Build AnnData (hd5) object from a counts table.",
        main=_build,
        options=[
            inputs, 
            sample_id,
            plot_format,
            outputs,
        ],
    )
    filter = CLICommand(
        "filter", 
        description="Filter AnnData (hd5) objects based on gene and cell counts.",
        main=_filter,
        options=[
            inputs, 
            min_counts_per_cell,
            min_cells_per_gene,
            min_cells,
            min_genes,
            plot_format,
            outputs,
        ],
    )
    cluster = CLICommand(
        "cluster", 
        description="Cluster cells based on their gene counts.",
        main=_cluster,
        options=[
            inputs, 
            plot_format,
            outputs,
        ],
    )
    
    app = CLIApp(
        'anndata-utils', 
        description='Analyze AnnData objects with scanpy.',
        version=__version__,
        commands=[
            build, 
            filter, 
            cluster,
        ],
    )

    app.run()

    return None


if __name__ == '__main__':
    main()