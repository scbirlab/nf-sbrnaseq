/*
 * Index the reference genome for use by STAR.
 */
process STAR_index {

   tag "${genome_id}"
   
   input:
   tuple val( genome_id ), path( fasta ), path( gff )

   output:
   tuple val( genome_id ), path( "star-index" )

   script:
   """
   mkdir star-index
   STAR \
      --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir star-index \
      --genomeFastaFiles ${fasta} \
      --sjdbOverhang 100 \
      --sjdbGTFfile ${gff} \
      --sjdbGTFfeatureExon CDS \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbGTFtagExonParentGene Parent \
      --sjdbGTFtagExonParentGeneName gene \
      --sjdbGTFtagExonParentGeneType gene_biotype

   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process STAR_align {

   tag "${id}:${genome_acc}" 
   label "big_cpu"
   time "2d"

   // errorStrategy 'retry'
   // maxRetries 1

   publishDir( 
      "${params.outputs}/mapped", 
      mode: 'copy',
      // saveAs: { "${id}.star.bam" },
      pattern: "*.bam"
   )
   publishDir( 
      "${params.outputs}/mapped", 
      mode: 'copy',
      // saveAs: { "${id}.${it}" },
      pattern: "*.{bg,tab}"
   )
   publishDir( 
      "${params.outputs}/mapped", 
      mode: 'copy',
      // saveAs: { "${id}.star.log" },
      pattern: "Log.final.out"
   )

   input:
   tuple val( id ), path( reads ), path( idx ), val( genome_acc )
   val strand

   output:
   tuple val( id ), path( "${id}.star.bam" ), emit: main
   tuple val( id ), path( "*.star.bg" ), emit: bg
   tuple val( id ), path( "*.star.tab" ), emit: tab
   path "${id}.star.log", emit: logs

   script:
   """
   STAR \
      --runMode alignReads \
      --runThreadN ${task.cpus} \
      --genomeDir "${idx}" \
      --readFilesIn ${reads} \
      --readFilesCommand zcat \
      --alignEndsType Local \
      --alignIntronMax 1 \
      --outFilterMultimapNmax 20 \
      --outSAMprimaryFlag AllBestScore \
      --outSAMattributes NH HI NM MD AS nM \
      --outSAMattrIHstart 0 \
      --twopassMode None \
      --quantMode GeneCounts \
      --outWigType bedGraph read1_5p \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingThreadN ${task.cpus} \
      --limitBAMsortRAM ${Math.round(task.memory.getBytes() * 0.8)}

   mv "Aligned.sortedByCoord.out.bam" "${id}.star.bam"
   mv "Log.final.out" "${id}.star.log"
   mv ReadsPerGene.out.tab "${id}.star.tab"

   for f in *.bg
   do
      mv "\$f" "\$(basename \$f .bg)".star.bg
   done
   
   """
}
