process join_featurecounts_UMItools {

   tag "${id}"
   label 'big_mem'
   cpus 1

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( umitools_table ), path( featurecounts_table )

   output:
   tuple val( id ), path( "all-counts.tsv" )

   script:
   """
   #!/usr/bin/env python

   import pandas as pd

   featurecounts_df = (
      pd.read_csv("${featurecounts_table}", sep="\\t", comment='#')
      .rename(columns={"Chr": "chr"})
      .assign(
         sample_id="${id}",
         gene_id=lambda x: x["Geneid"],
      )
   )

   bam_cols = [col for col in featurecounts_df if col.endswith(".bam")]
   featurecounts_df = featurecounts_df.rename(
      columns={col: "pseudobulk_read_count" for col in bam_cols}
   )

   df_in = (
      pd.read_csv("${umitools_table}", sep="\\t")
      .assign(
         sample_id="${id}",
      )
   )

   df_out = (
      featurecounts_df
      .merge(
         df_in,
         how="outer",
      )
      .drop_duplicates()
      .assign(
         umi_count=lambda x: x["umi_count"].fillna(0).astype(int),
         read_count=lambda x: x["read_count"].fillna(0).astype(int),
      )
   )

   #if df_out.shape[0] != df_in.shape[0]:
   #   raise ValueError(
   #      f"Extra rows were added! {df_in.shape[0]=} -> {df_out.shape[0]=} "
   #      f"(Featurecounts had {featurecounts_df.shape[0]} rows)."
   #   )

   df_out.to_csv("all-counts.tsv", sep="\\t", index=False)

   """
}

process count_genomes_per_cell {

   tag "${id}"
   label "big_mem"

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
      // pattern: "genome-per-cell*.tsv",
   )

   input:
   tuple val( id ), path( joined_table )

   output:
   tuple val( id ), path( "genome-per-cell{,-summary}.tsv" ), emit: chr_per_cell
   tuple val( id ), path( "*.{png,tsv}" ), emit: plots, optional: true

   script:
   """
   #!/usr/bin/env python

   import pandas as pd
   from carabiner.mpl import figsaver, scattergrid

   GENOME_KEY = "genome_accession"
   CELL_BC_KEY = "cell_barcode"

   df = pd.read_csv("${joined_table}", sep="\\t")#.query("not gene_biotype.isin(('rRNA', 'tRNA'))")
   n_genomes = df[GENOME_KEY].nunique()

   df_sum = (
      df
      .groupby([CELL_BC_KEY, GENOME_KEY])
      [["umi_count", "read_count"]]
      .sum()
      .reset_index()
   )

   cell_sums = (
      df
      .groupby(CELL_BC_KEY)
      [["umi_count", "read_count"]]
      .sum()
   )

   m = (
      df_sum
      .pivot(
         index=CELL_BC_KEY,
         columns=GENOME_KEY,
         values=["umi_count", "read_count"],
      )
      .fillna(0)
      .astype(int)
   )
   print(m.head())
   m = m.reindex(cell_sums.query("read_count >= 2").index)
   m.columns = [":".join(c) for c in m.columns.to_flat_index()]
   m.to_csv("genome-per-cell.tsv", sep="\\t")

   if len(m.shape) > 1 and m.shape[-1] > 0:
      fig, axes = scattergrid(
         m,
         grid_columns=m.columns.tolist(),
         #log=m.columns.tolist(),
      )
      figsaver(format="png")(
         fig=fig,
         name="genome-per-cell",
         df=m.reset_index(),
      )

      df_hist = (
         df_sum
         .query("read_count >= 2")
         .groupby(CELL_BC_KEY)
         [[GENOME_KEY]]
         .nunique()
         .reset_index()
         .groupby(GENOME_KEY)
         [[CELL_BC_KEY]]
         .nunique()
         .reset_index()
         .rename(columns={
            GENOME_KEY: "n_genomes_per_cell",
            CELL_BC_KEY: "n_cells_with_n_genomes",
         })
      )

      df_hist.to_csv("genome-per-cell-summary.tsv", sep="\\t", index=False)

   """

}
