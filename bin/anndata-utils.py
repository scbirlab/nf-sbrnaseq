#!/usr/bin/env python

"""Ultilities for building AnnData objects and basic clustering analysis."""

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

from scipy.sparse import csr_matrix, coo_matrix, issparse
from pandas.api.types import CategoricalDtype

CELL_ID = ["cell_barcode"]
GENE_ID = ["genome_accession", "assembly", "chr", "gene_id"]
UMI_COUNT_COLUMN = "umi_count"
MIN_COUNTS_PER_CELL = 60
MIN_CELLS_PER_GENE = 3
MIN_GENES_TO_KEEP = 200
MIN_CELLS_TO_KEEP = 950

__version__ = "0.0.1"

def _joiner(col: str, char: str = ":"):
    def f(df: pd.DataFrame) -> pd.Series:
        return df[col].astype(str).apply(
            char.join, 
            raw=True, 
            axis=1, 
            result_type="reduce",
        )
    return f


def make_anndata_from_table(counts_table: str) -> pd.DataFrame:

    df = (
        pd.read_csv(
            counts_table, 
            sep="\t", 
            low_memory=False,
        )
        .query(f"{UMI_COUNT_COLUMN} > 0")
        .assign(**{
            "__cell_id__": _joiner(CELL_ID),
            "__gene_id__": _joiner(GENE_ID),
        })
    )
    
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

    matrix_df = pd.DataFrame.sparse.from_spmatrix(
        matrix,
        index=cell_cat.categories,
        columns=gene_cat.categories,
    )
    adata = ad.AnnData(matrix_df)
    adata.obs_names = matrix_df.index
    adata.var_names = matrix_df.columns

    print_err("shape:", matrix_df.shape)
    print_err("dup genes?:", matrix_df.index.duplicated().any())
    print_err("n dup genes:", matrix_df.index.duplicated().sum())

    print_err("dup cells?:", matrix_df.columns.duplicated().any())
    print_err("n dup cells:", matrix_df.columns.duplicated().sum())

    print_err("var has duplicates?:", adata.var_names.duplicated().any())
    print_err("obs has duplicates?:", adata.obs_names.duplicated().any())

    return df, adata


def save_anndata(x, filename: str, printf: Callable = print) -> None:
    printf(f"Saving the following AnnData object as {filename}.h5ad:")
    printf(x)
    x.write(f"{filename}.h5ad", compression="gzip")
    return None


def plot_barnyard(
    adata,
    filename: str,
    grouping: str = "genome_accession",
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
    count_df, adata = make_anndata_from_table(args.input)
    adata.obs["sample_id"] = args.id

    print_err(count_df)
    print_err(adata)

    print_err("Annotating gene features...")
    gene_ann_columns = GENE_ID + [
        "locus_tag", 
        "Name", 
        "gene_biotype",
        "Length",
        "pseudobulk_read_count",
    ]
    gene_ann_df = (
        count_df[["__gene_id__"] + gene_ann_columns]
        .drop_duplicates()
        .set_index("__gene_id__")
    )
    print_err(gene_ann_df.shape, sorted(gene_ann_df.columns))
    gene_ann_df = gene_ann_df.reindex(adata.var.index)
    print_err(gene_ann_df)
    print_err(adata.var)
    print_err(sorted(set(adata.var.index) - set(gene_ann_df.index)))
    d = {idx: len([idx2 for idx2 in gene_ann_df.index if idx2==idx]) for idx in gene_ann_df.index}
    print_err(d)
    print_err(sorted(d, key=lambda x: d[x]))
    adata.var = adata.var.join(gene_ann_df, how="left")
    # for col in gene_ann_columns:
    #     adata.var[col.casefold()] = gene_ann_df[col]

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
            "Length", 
            "pct_dropout_by_counts",
            "pseudobulk_read_count",
        ],
        log=[
            "n_cells_by_counts", 
            "total_counts",
            "pseudobulk_read_count",
        ],
        grouping="genome_accession",
        aspect_ratio=1.25,
    )
    _figsaver(
        fig=fig,
        name=f"{args.output}.gene-qc",
        df=adata.var,
    )

    adata.layers["counts"] = adata.X.copy()

    save_anndata(adata, f"{args.output}", printf=print_err)

    return None

    
@clicommand(message="Filtering AnnData with the following parameters")
def _filter(args: Namespace) -> None:

    _figsaver = figsaver(format=args.plot_format)
    adata = sc.read_h5ad(args.input)

    if "counts" not in adata.layers:
        raise ValueError("Missing raw counts layer")

    if args.min_counts_per_cell < 1 or args.min_cells_per_gene < 1:
        raise ValueError("Filtering thresholds must be positive integers")

    n_cells_before, n_genes_before = adata.shape

    # Use untransformed molecule counts for filtering.
    adata.X = adata.layers["counts"].copy()

    cell_counts = np.asarray(
        adata.layers["counts"].sum(axis=1)
    ).ravel()
    # Retention curve independent of the chosen threshold.
    thresholds = np.unique(np.concatenate((
        np.array([0, args.min_counts_per_cell]),
        np.arange(0, max(101, int(cell_counts.max()) + 2), 10),
    )))

    retained = np.array([
        np.count_nonzero(cell_counts >= threshold)
        for threshold in thresholds
    ])

    retention = pd.DataFrame({
        "min_counts_per_cell": thresholds,
        "n_cells": retained,
        "fraction_retained": retained / n_cells_before,
    })
    
    fig, ax = grid(aspect_ratio=1.5)
    ax.plot(thresholds, retained)
    ax.axvline(
        args.min_counts_per_cell,
        color="lightgrey",
    )
    ax.set(
        xlabel="Minimum UMIs per cell",
        ylabel="Cells retained",
    )
    _figsaver(
        fig=fig,
        name=f"{args.output}.cell-retention",
        df=retention,
    )

    if adata.n_obs == 0:
        raise ValueError(
            "No cells passed the UMI threshold"
        )

    # Apply explicit cell threshold.
    keep_cells = cell_counts >= args.min_counts_per_cell
    adata = adata[keep_cells, :].copy()

    if adata.n_obs == 0:
        raise ValueError("No cells passed the UMI threshold")

    # Filter genes based on detection among retained cells.
    sc.pp.filter_genes(
        adata,
        min_cells=args.min_cells_per_gene,
    )
    # Remove cells left empty after gene filtering.
    remaining_counts = np.asarray(
        adata.layers["counts"].sum(axis=1)
    ).ravel()

    n_empty = int(np.count_nonzero(remaining_counts == 0))
    print_err(f"Cells emptied by gene filtering: {n_empty}")

    adata = adata[remaining_counts > 0].copy()

    if adata.n_vars == 0:
        raise ValueError(
            "No genes passed the detection threshold"
        )
    remaining = np.asarray(
        adata.layers["counts"].sum(axis=1)
    ).ravel()

    adata = adata[remaining > 0, :].copy()

    if adata.n_obs == 0:
        raise ValueError("No cells with nonzero counts after filtering")
    # Record the applied thresholds and retention.
    adata.uns["filtering"] = {
        "min_counts_per_cell": args.min_counts_per_cell,
        "min_cells_per_gene": args.min_cells_per_gene,
        "n_cells_before": n_cells_before,
        "n_cells_after": adata.n_obs,
        "n_genes_before": n_genes_before,
        "n_genes_after": adata.n_vars,
    }

    print_err(
        f"Retained {adata.n_obs}/{n_cells_before} cells "
        f"and {adata.n_vars}/{n_genes_before} genes"
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

    # Preserve original molecule counts for downstream analyses.
    if "counts" not in adata.layers:
        raise ValueError("Missing raw counts layer")

    # Ensure only cells with positive total counts are normalised.
    library_sizes = np.asarray(
        adata.layers["counts"].sum(axis=1)
    ).ravel()
    if np.any(library_sizes <= 0):
        raise ValueError(
            "Cannot normalise cells with zero molecule counts"
        )

    # Standard library-size normalisation for exploratory
    # dimensionality reduction.
    adata.var["not_is_rRNA"] = (~adata.var["is_rRNA"])
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    adata.uns["normalization"] = {
        "method": "normalize_total_log1p",
        "raw_counts_layer": "counts",
    }

    print_err("Annotating highly variable genes...")
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
        use_highly_variable = False #adata.var['highly_variable_non_rRNA'].sum() >= 3

    
    print_err("Applying PCA...")

    mask = (
        "highly_variable_non_rRNA"
        if use_highly_variable else "not_is_rRNA"
    )

    if int(adata.var[mask].sum()) < 3:
        raise ValueError(
            "Fewer than three non-rRNA genes available for PCA"
        )

    n_comps = min(50, adata.n_obs - 1, int(adata.var[mask].sum()) - 1)

    if n_comps < 2:
        raise ValueError("Insufficient dimensions for PCA")

    sc.pp.pca(
        adata,
        n_comps=n_comps,
        mask_var=mask,
        random_state=42,
    )

    print_err("Getting nearest neighbours...")
    sc.pp.neighbors(
        adata,
        n_neighbors=min(15, adata.n_obs - 1),
        metric="euclidean",
        use_rep="X_pca",
        random_state=42,
    )

    print_err("Carrying out Leiden clustering...")
    sc.tl.leiden(
        adata,
        key_added="leiden",
        resolution=1.0,
        flavor="igraph",
        directed=False,
        random_state=42,
    )

    print_err(
        f"Found {adata.obs['leiden'].cat.categories.size} clusters"
    )

    if adata.obs["leiden"].nunique() > 1:
        sc.tl.rank_genes_groups(
            adata,
            "leiden",
            mask_var=mask,
            method="t-test",
            rankby_abs=True,
            use_raw=False,
        )

    print_err("Embedding as UMAP...")
    sc.tl.umap(
        adata,
        min_dist=.1,
        spread=1.,
        random_state=42,
    )

    save_anndata(adata, f"{args.output}", printf=print_err)

    n_genes, n_cells = adata.n_vars, adata.n_obs
    biotype_flags = [c for c in adata.obs if c.startswith("is_")]
    colors_to_plot = [
        "total_counts", 
        "n_genes_by_counts",
        "doublet_score",
        "predicted_doublet",
        "leiden",
        "tRNA_rRNA_ratio",
        "pseudobulk_read_count",
    ] + [
        c for c in adata.obs if c.startswith("pct_counts_")
    ]
    colors_to_plot = [c for c in colors_to_plot if c in adata.obs]
    log_colors = [c for c in colors_to_plot if c.endswith("_counts")] + [
        "tRNA_rRNA_ratio", 
        "pseudobulk_read_count",
    ]

    print_err("Making UMAP plots...")
    fig, axes = grid(
        ncol=len(colors_to_plot), 
        aspect_ratio=1., 
        panel_size=3.75,
    
    )
    for ax, col in zip(axes, colors_to_plot):
        colors = adata.obs[col]
        col_is_categorical = isinstance(colors, pd.Series) and colors.dtype in ("category", "object", "str")
        use_log = (not col_is_categorical and np.all(colors > 0.) and col in log_colors)
        plotter = partial(
            ax.scatter,
            s=.5,
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
            bbox_to_anchor=(.5, -.2),
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
        df=pd.DataFrame(
            adata.obsm['X_umap'], 
            columns=["UMAP_1", "UMAP_2"], 
            index=adata.obs_names,
        ),
    )
    print_err("Plotting heatmap of mean gene markers per cell cluster...")
    # Compute a dendrogram on a dense representation only (avoids pandas+csr bug)
    try:
        if "X_pca" in adata.obsm:
            sc.tl.dendrogram(adata, groupby="leiden", use_rep="X_pca", key_added="dendrogram_leiden")
            dendro_arg = "dendrogram_leiden"
        else:
            if issparse(adata.X):
                adata.obsm["X_dense"] = adata.X.toarray()
                rep = "X_dense"
            else:
                rep = "X"
            sc.tl.dendrogram(adata, groupby="leiden", use_rep=rep, key_added="dendrogram_leiden")
            dendro_arg = "dendrogram_leiden"
    except Exception as e:
        print(f"[warn] dendrogram computation failed: {e}. Proceeding without dendrogram.")
        dendro_arg = False

    if adata.obs["leiden"].nunique() > 1:
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