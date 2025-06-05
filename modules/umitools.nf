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


   container params.umitools_image

   input:
   tuple val( id ), path( reads ), val( umis ), path( whitelist )

   output:
   tuple val( id ), path( "*.extracted.fastq.gz" ), emit: main
   path "*.log", emit: logs

   script:
   def error_correct = ( params.allow_cell_errors ? "--error-correct-cell" : "" )
   def second_reads  = ( 
      reads[1] 
      ? "--bc-pattern '${umis[0]}' --bc-pattern2 '${umis[1]}' --read2-in '${reads[1]}' --read2-out \${f}_R2.fastq.gz"
      : "--bc-pattern '${umis}'" 
   )
   """
   MAX_SIZE=500000000  # 500 million
   split -l \$MAX_SIZE "${whitelist}" chunk_

   for f in chunk_*
   do
      umi_tools extract \
         --extract-method regex ${error_correct} \
         --quality-encoding phred33 \
         --whitelist "\$f" \
         --log "${id}.\$f.extract.log" \
         --stdin "${reads[0]}" \
         --stdout "\$f"_R1.fastq.gz ${second_reads}
   done

   for i in 1 ${reads[1] ? '2' : ''}
   do
      cat chunk_*_R\$i.fastq.gz > ${id}_R\$i.extracted0.fastq.gz
   done

   # make sure no 0-length reads
   for f in ${id}_R?.extracted0.fastq.gz
   do
      zcat \$f \
         | awk \
            '
            NR%4 == 1 { name = \$0 } 
            NR%4 == 2 { if (length(\$0) > 0) { seq = \$0 } else { seq = "N" } } 
            NR%4 == 0 { if (length(\$0) > 0) { qual = \$0 } else { qual = "?" } } 
            ( NR > 1 && NR%4 == 1 ) { print name ORS seq ORS "+" ORS qual } 
            END { print name ORS seq ORS "+" ORS qual }
            ' \
         | gzip \
      > \$(basename \$f .fastq.gz).extracted.fastq.gz
   done
   rm ${id}_R?.extracted0.fastq.gz chunk_*_R?.fastq.gz

   """
}

/*
 * Count unique UMIs per cell per gene
 */
process UMItools_count {

   tag "${sample_id}"

   // errorStrategy 'retry'
   // maxRetries 2

   container params.umitools_image

   input:
   tuple val( sample_id ), path( bamfile )

   output:
   tuple val( sample_id ), path( "*.tsv" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   samtools sort ${bamfile} -o ${bamfile.getBaseName()}.sorted.bam
   samtools index ${bamfile.getBaseName()}.sorted.bam
   umi_tools count \
      --paired \
		--per-gene \
      --per-cell \
		--gene-tag XT \
      --chimeric-pairs discard \
      --unpaired-reads discard \
      --log ${sample_id}.count.log \
		--stdin ${bamfile.getBaseName()}.sorted.bam \
      --stdout ${sample_id}.umitools_count.tsv \
   && rm ${bamfile.getBaseName()}.sorted.bam

   """
}