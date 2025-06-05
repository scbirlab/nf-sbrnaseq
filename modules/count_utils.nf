process join_featurecounts_UMItools {

   tag "${id}"
   label 'some_mem'

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
   get_col_number () (
      local col="\$1"
      local sep=\${2:-,}
      head -n1 | tr "\$sep" \$'\\n' | grep -n "\$col" | cut -d: -f1
   )
   sort_table () (
      local col=\${1:-"1"}
      cat > "temp-table"
      head -n1 "temp-table" \
      | cat - <(tail -n+2 "temp-table" \
      | LC_ALL=C sort -k"\$col" -d -t, --parallel=${task.cpus} --buffer-size="${task.memory.getBytes()}b") && rm "temp-table"
   )
   set -x
   grep -v '^#' "${featurecounts_table}" | tr \$'\\t' , > "featurecounts-table0.csv"
   LOCUS_TAG_COL=\$(get_col_number "Geneid" < "featurecounts-table0.csv")
   # rename header, add locus tag column, and sort
   sed 's/,${id}.*\\.umicollapse\\.bam\$/,bulk_read_count/' \
      "featurecounts-table0.csv" \
   | awk -F, -v OFS=, \
      -v locus_col="\$LOCUS_TAG_COL" -v sample_id="${id}" \
      'NR == 1 { print "gene_id", "sample_id", \$0 } NR > 1 { print \$locus_col, sample_id, \$0 } ' \
   | sort_table \
   > "featurecounts-table.csv"

   cat ${umitools_table} | tr \$'\\t' , > "umitools-table0.csv"
   sed 's/,count\$/,umi_count/;s/^gene,/gene_id,/;s/,cell,/,cell_barcode,/' \
      "umitools-table0.csv" \
   | sort_table \
   > "umitools-table.csv"

   LC_COLLATE=C join --header -j1 -t, --nocheck-order \
      "featurecounts-table.csv" "umitools-table.csv" \
   | tr , \$'\t' \
   > "all-counts.tsv" 
   """
}

process count_genomes_per_cell {

   tag "${id}"

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
      pattern: "genome-per-cell*.tsv",
   )

   input:
   tuple val( id ), path( joined_table )

   output:
   tuple val( id ), path( "genome-per-cell{,-summary}.tsv" ), emit: chr_per_cell

   script:
   """
   get_col_number () (
      local col="\$1"
      local sep=\${2:-,}
      head -n1 | tr "\$sep" \$'\\n' | grep -n "\$col" | cut -d: -f1
   )
   set -x
   
   # count the number of assemblies per cell
   cell_bc_col=\$(get_col_number "cell_barcode" \$'\t' < "${joined_table}")
   chr_col=\$(get_col_number "Chr"  \$'\t' < "${joined_table}")
   genome_col=\$(get_col_number "genome_accession"  \$'\t' < "${joined_table}")
   awk -F'\t' -v OFS='\t' \
      -v bc_col="cell_barcode" -v g_col="genome_accession" \
      '
      BEGIN { print "cell_barcode", "n_genomes" } 
      NR == 1 { for (i = 1; i <= NF; i++) col_idx[\$i]=i } 
      NR > 1 { a[\$col_idx[bc_col]][\$col_idx[g_col]]++; if (a[\$col_idx[bc_col]][\$col_idx[g_col]] == 1) counts[\$col_idx[bc_col]]++ } 
      END { for (bc in counts) print bc, counts[bc] }
      ' \
      "${joined_table}" \
   > "genome-per-cell.tsv"

   awk -F'\t' -v OFS='\t' \
      '
      BEGIN { print "n_genomes", "n_cell_barcodes_with_n_genomes" } 
      NR > 1 { a[\$2]++ } 
      END { for (key in a) print key, a[key] }
      ' \
      "genome-per-cell.tsv" \
   > "genome-per-cell-summary.tsv"
   """
}