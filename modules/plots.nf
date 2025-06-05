process plot_UMI_distributions {

   tag "${id}"
   label 'big_mem'
   cpus 1

   publishDir( 
      "${params.outputs}/plot-counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( umi_table )

   output:
   tuple val( id ), path( "*.{png,csv}" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   biotype_blocklist = [
      #"SRP_RNA",
      #"antisense_RNA",
   ]
   df = (
      pd.read_csv(
         "${umi_table}", 
         sep='\\t',
         low_memory=False,
      )
      .query("~gene_biotype.isin(@biotype_blocklist)")  # very low representation, breaks histogram bins
   )

   fig, axes = scattergrid(
      df,
      grid_columns=["umi_count", "bulk_read_count"],
      log=["umi_count", "bulk_read_count"],
      grouping=["gene_biotype"],
      aspect_ratio=1.25,
   )
   #for ax in fig.axes:
   #   ax.set(yscale="log")
   figsaver(format="png")(
      fig=fig,
      name='umi-hist',
      df=df,
   )

   """
   
}


process plot_cells_per_gene_distribution {

   tag "${id}"
   label 'big_mem'
   cpus 1

   publishDir( 
      "${params.outputs}/plot-counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( umi_table )

   output:
   tuple val( id ), path( "*.{png,csv}" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   df = (
      pd.read_csv(
         "${umi_table}", 
         sep='\\t',
      ) 
      .groupby(["genome_accession", "gene_id"])
   )
   bc_counts = (
      df[["cell_barcode"]]
      .count()
      .reset_index()
      .merge(
         df[["umi_count", "bulk_read_count"]]
         .sum()
         .reset_index()
      )
      .rename(columns={"cell_barcode": "cell_count"})
   )

   fig, axes = scattergrid(
      bc_counts,
      grid_columns=["umi_count", "bulk_read_count", "cell_count"],
      log=["umi_count", "bulk_read_count", "cell_count"],
      grouping=["genome_accession"],
      aspect_ratio=1.25,
   )
   figsaver(format="png")(
      fig=fig,
      name='umis-and-cells-per-gene',
      df=bc_counts,
   )

   """
   
}


process plot_genes_per_cell_distribution {

   tag "${id}"
   label 'big_mem'
   cpus 1

   publishDir( 
      "${params.outputs}/plot-counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( umi_table )

   output:
   tuple val( id ), path( "*.{png,csv}" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd

   df = pd.read_csv(
      "${umi_table}", 
      sep='\\t',
   ) 
   renamer = {"gene_id": "gene_count"}

   genes_per_cell = (
      df
      .groupby("cell_barcode")[["gene_id"]]
      .count()
      .reset_index()
      .merge(
         df
         .groupby("cell_barcode")[["umi_count"]]
         .sum()
         .reset_index()
      )
      .rename(columns=renamer)
      .sort_values("umi_count")
   )

   genes_per_biotype_per_cell = (
      df
      .groupby(["cell_barcode", "gene_biotype"])[["gene_id"]]
      .count()
      .reset_index()
      .merge(
         df
         .groupby(["cell_barcode", "gene_biotype"])[["umi_count"]]
         .sum()
         .reset_index()
      )
      .rename(columns=renamer)
      .sort_values("umi_count")
   )

   
   for grouping in (
      ["cell_barcode"],
      ["cell_barcode", "gene_biotype"],
   ):
      grouped_df = df.groupby(grouping)
      count_df = (
         grouped_df[["gene_id"]]
         .count()
         .reset_index()
         .merge(
            grouped_df[["umi_count", "bulk_read_count"]]
            .sum()
            .reset_index()
         )
         .rename(columns=renamer)
         .sort_values("umi_count")
      )
      using_biotype = "gene_biotype" in grouping
      fig, axes = scattergrid(
         count_df,
         grid_columns=["umi_count", "bulk_read_count", "gene_count"],
         log=["umi_count", "bulk_read_count", "gene_count"],
         grouping=["gene_biotype"] if using_biotype else None,
         aspect_ratio=1.25 if using_biotype else 1.,
      )
      figsaver(format="png")(
         fig=fig,
         name="-".join(grouping),
         df=count_df,
      )

   """
   
}