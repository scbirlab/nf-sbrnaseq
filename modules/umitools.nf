/*
 * Extract cell barcodes and UMIs
 */
process UMItools_extract {

   tag "${id}"
   label "big_mem"
   time "2d"  // in case of very large whitelists

   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/extracted", 
      mode: 'copy',
      pattern: "*.extract.log", 
   )

   input:
   tuple val( id ), path( reads ), val( umis ), path( whitelist )

   output:
   tuple val( id ), path( "*.extracted.fastq.gz" ), emit: main
   path "*.log", emit: logs

   script:
   def error_correct = ( params.allow_cell_errors ? "--error-correct-cell" : "" )
   def second_reads  = ( 
      reads[1]
      ? "--bc-pattern2 '${umis[1]}' --read2-in '${reads[1]}' --read2-out '${id}'_R2.extracted0.fastq.gz"
      : "" 
   )
   """
   zcat "${whitelist}" > wl.txt
   umi_tools extract \
      --extract-method regex ${error_correct} \
      --quality-encoding phred33 \
      --whitelist "wl.txt" \
      --bc-pattern '${umis[0]}' \
      --log "${id}.${whitelist.baseName}.extract.log" \
      --stdin "${reads[0]}" \
      --stdout "${id}"_R1.extracted0.fastq.gz ${second_reads}

   # make sure no 0-length reads
   for f in "${id}"_R?.extracted0.fastq.gz
   do
      zcat \$f \
      | awk '
         NR % 4 == 1 { name = \$0 } 
         NR % 4 == 2 { if (length(\$0) > 0) { seq = \$0 } else { seq = "N" } } 
         NR % 4 == 0 { if (length(\$0) > 0) { qual = \$0 } else { qual = "?" } } 
         ( NR > 1 && NR % 4 == 1 ) { print name ORS seq ORS "+" ORS qual } 
         END { print name ORS seq ORS "+" ORS qual }
         ' \
      | pigz -v --stdout -p ${task.cpus} \
      > \$(basename \$f .fastq.gz).extracted.fastq.gz
   done

   rm "${id}"_R?.extracted0.fastq.gz

   """
}

process concat_extractions {

   tag "${id}"

   // errorStrategy 'retry'
   // maxRetries 2

   // publishDir( 
   //    "${params.outputs}/extracted", 
   //    mode: 'copy',
   //    pattern: "*.extract.log", 
   // )

   input:
   tuple val( id ), path( reads, stageAs: '??/*' )

   output:
   tuple val( id ), path( "*.extracted-concat.fastq.gz", arity: 1..2 )

   script:
   """
   cat ??/*_R1.*extracted.fastq.gz > "${id}"_R1.extracted-concat.fastq.gz
   if ls ??/*_R2.*extracted.fastq.gz 1> /dev/null 2>&1
   then
      cat ??/*_R2.*extracted.fastq.gz > "${id}"_R2.extracted-concat.fastq.gz
   fi

   """

}

process bam2table {

   tag "${id}"

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( bamfile )
   val paired

   output:
   tuple val( id ), path( "tab.tsv" )

   script:
   """
   samtools view -h -F${paired ? "1664 -f3" : "1536"} \
      --tag 'XS:Assigned' \
      "${bamfile}" \
      -bS -o filtered.bam

   samtools view filtered.bam \
   | awk -F'\\t' -v OFS='\\t' '
      BEGIN { print "read_id", "gene_id", "chr", "cell_barcode", "umi", "umi_cluster_size", "umi_read_count" }
      {
         assign=""; gene=""; cb=""; umi=""; umi_read_count=""; umi_cluster_size="";
         chr=\$3;
         delete b; delete a;
         split(\$1, b, "_");
         cb = b[length(b) - 1]
         for (i=12; i<=NF; i++) {
            if (\$i ~ /^XS:Z:/) { split(\$i, a, ":"); assign = a[3] }
            else if (\$i ~ /^XT:Z:/) { split(\$i, a, ":"); gene = a[3] }
            else if (\$i ~ /^RX:Z:/) { split(\$i, a, ":"); umi = a[3] }
            else if (\$i ~ /^cs:i:/) { split(\$i, a, ":"); umi_cluster_size = a[3] }
            else if (\$i ~ /^su:i:/) { split(\$i, a, ":"); umi_read_count = a[3] }
         }
         if ( assign != "Assigned" || gene == "" || umi == "" || umi_read_count == "" ) next;
         print \$1, gene, chr, cb, umi, umi_cluster_size, umi_read_count
      }
   ' \
   > tab.tsv

   rm filtered.bam

   """
}


process counts_per_gene_awk {

   tag "${id}"

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "counts_per_gene.tsv" )

   script:
   """
   awk -F'\\t' -v OFS='\\t' '
      BEGIN { print "chr", "gene_id", "umi_count", "read_count" } 
      NR > 1 {
         u1[\$3,\$2]++; 
         r[\$3,\$2] += \$7 
      } 
      END { 
         for (key in u) {
            split( key, k, SUBSEP ); 
            print k[1], k[2], u[key], r[key]
         }
      }
   '  "${tabfile}" \
   > "counts_per_gene.tsv"

   """
}


process counts_per_gene {

   tag "${id}"
   label "big_mem"
   cpus 1

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "counts_per_gene.tsv" )

   script:
   """
   #!/usr/bin/env python
   import pandas as pd

   (
      pd.read_csv(
         "${tabfile}",
         sep="\\t",
      )
      .groupby(["chr", "gene_id"])
      .agg({
         "umi": "nunique",
         "umi_read_count": "sum",
      })
      .rename(columns={
         "umi": "umi_count",
         "umi_read_count": "read_count",
      })
      .rename(columns={
         "umi": "umi_count",
      })
      .to_csv(
         "counts_per_gene.tsv",
         sep="\\t",
         index=True,
      )
   )
   
   """
}



process counts_per_cell_awk {

   tag "${id}"

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "counts_per_cell.tsv" )

   script:
   """
   awk -F'\\t' -v OFS='\\t' '
      BEGIN { print "cell_barcode", "umi_count", "read_count" } 
      NR>1 { u[\$4]++; r[\$4] += \$7 } 
      END { 
         for (key in u) { 
            print key, u[key], r[key]
         }
      }
   '  "${tabfile}" \
   > "counts_per_cell.tsv"

   """
}


process counts_per_cell {

   tag "${id}"
   label "big_mem"
   cpus 1

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "counts_per_cell.tsv" )

   script:
   """
   #!/usr/bin/env python
   import pandas as pd

   (
      pd.read_csv(
         "${tabfile}",
         sep="\\t",
      )
      .groupby("cell_barcode")
      .agg({
         "umi": "nunique",
         "umi_read_count": "sum",
      })
      .rename(columns={
         "umi": "umi_count",
         "umi_read_count": "read_count",
      })
      .to_csv(
         "counts_per_cell.tsv",
         sep="\\t",
         index=True,
      )
   )

   """
}


process counts_per_gene_per_cell_awk {

   tag "${id}"
   label "big_mem"
   cpus 1

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "counts_per_gene_per_cell.tsv" )

   script:
   """
   awk -F'\\t' -v OFS='\\t' '
      BEGIN { print "chr", "gene_id", "cell_barcode", "umi_count", "read_count" } 
      NR > 1 && NF >= 7 && \$7 > 0 { u[\$3,\$2,\$4]++; r[\$3,\$2,\$4] += \$7 } 
      END { 
         for (key in u) {
            if ( r[key] > 0) {
               split( key, k, SUBSEP ); 
               print k[1], k[2], k[3], u[key], r[key]
            }
         }
      }
   '  "${tabfile}" \
   > "counts_per_gene_per_cell.tsv"

   """
}


process counts_per_gene_per_cell {

   tag "${id}"
   label "big_mem"
   cpus 1

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "counts_per_gene_per_cell.tsv" )

   script:
   """
   #!/usr/bin/env python
   import pandas as pd

   (
      pd.read_csv(
         "${tabfile}",
         sep="\\t",
      )
      .groupby(["chr", "gene_id", "cell_barcode"])
      .agg({
         "umi": "nunique",
         "umi_read_count": "sum",
      })
      .rename(columns={
         "umi": "umi_count",
         "umi_read_count": "read_count",
      })
      .to_csv(
         "counts_per_gene_per_cell.tsv",
         sep="\\t",
         index=True,
      )
   )

   """
}


// Count unique UMIs per cell per gene
process UMItools_count_tab {

   tag "${id}"

   label 'big_mem'

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}-${it}" },
   )

   input:
   tuple val( id ), path( tabfile )
   val umi_method
   val per_clone

   output:
   tuple val( id ), path( "${id}.umitools_count.tsv" ), emit: main
   path "${id}.umitools_count.log", emit: logs

   script:
   """
   export MPLCONFIGDIR=tmp
   mkdir \$MPLCONFIGDIR

   umi_tools count_tab \
      --method ${umi_method} ${per_clone ? "--per-cell" : ""} \
		--stdin "${tabfile}" \
      --stdout umitools_count0.tsv \
      --log "${id}.umitools_count.log"

   awk -v OFS='\\t' -v id="${id}" '
      BEGIN { print "sample_id", "guide_name", "umi_count" }
      NR > 1 { \$1 = \$1; print id, \$0 }
   ' umitools_count0.tsv \
   > "${id}.umitools_count.tsv"

   """
}


/*
 * Count unique UMIs per cell per gene
 */
process UMItools_count {

   tag "${id}"

   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      // saveAs: { "${id}.${it}" }, 
   )

   input:
   tuple val( id ), path( bamfile )
   val paired

   output:
   tuple val( id ), path( "*.umitools_count.tsv" ), emit: main
   path "*.umitools_count.log", emit: logs

   script:
   """
   samtools view -F1024 "${bamfile}" \
   | samtools sort - -@ ${task.cpus} -m ${Math.round(task.memory.getGiga() * 0.8)}G -o sorted.bam
   samtools index sorted.bam
   umi_tools count \
		--per-gene \
      --per-cell \
		--gene-tag XT \
      --assigned-status-tag XS ${paired ? '--paired --chimeric-pairs discard --unpaired-reads discard' : ''} \
      --log ${id}.umitools_count.log \
		--stdin sorted.bam \
      --stdout ${id}.umitools_count.tsv \
   && rm sorted.bam

   """
}
