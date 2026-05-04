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

   from carabiner import print_err
   import pandas as pd

   featurecounts_df = (
      pd.read_csv(
         "${featurecounts_table}", 
         sep="\\t", 
         comment="#",
      )
      .rename(columns={
         "Chr": "chr", 
         "Geneid": "gene_id",
      })
      .assign(
         sample_id="${id}",
      )
   )

   bam_cols = [
      col for col in featurecounts_df 
      if col.endswith(".bam")
   ]
   if len(bam_cols) > 0:
      bam_cols = bam_cols[0]
   else:
      raise AttributeError(
         f"No column ending in .bam in featurecounts table: {featurecounts_df.columns}"
      )
   featurecounts_df = featurecounts_df.rename(
      columns={bam_cols: "pseudobulk_read_count"},
   )

   df_in = (
      pd.read_csv("${umitools_table}", sep="\\t")
      .assign(
         sample_id="${id}",
      )
   )

   print_err(
      f"Merging on {set(featurecounts_df.columns).intersection(set(df_in.columns))}"
   )
   df_out = (
      featurecounts_df
      .merge(
         df_in,
         on=["gene_id", "chr"],
         how="right",
         validate="many_to_one",
      )
      .drop_duplicates()
      .assign(
         umi_count=lambda x: x["umi_count"].fillna(0).astype(int),
         read_count=lambda x: x["read_count"].fillna(0).astype(int),
      )
   )

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
