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
      saveAs: { "${id}.${it}" },
      pattern: "*.bam"
   )
   publishDir( 
      "${params.outputs}/mapped", 
      mode: 'copy',
      saveAs: { "${id}.star.log" },
      pattern: "Log.final.out"
   )

   input:
   tuple val( id ), path( reads ), path( idx ), val( genome_acc )

   output:
   tuple val( id ), path( "Aligned.out.bam" ), emit: main
   path "Log.final.out", emit: logs

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
      --outFilterMultimapNmax 10 \
      --outSAMmapqUnique 255 \
      --outSAMprimaryFlag AllBestScore \
      --outSAMattributes All \
      --outSAMattrIHstart 0 \
      --twopassMode None \
      --quantMode GeneCounts \
      --outSAMtype BAM Unsorted

   """
}
