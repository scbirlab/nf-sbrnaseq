/*
 * Index the reference genome for use by Bowtie2.
 */
process bowtie2_index {

   tag "${id}"
   
   input:
   tuple val( id ), path( fasta )

   output:
   tuple val( id ), path( "*.bt2" )

   script:
   """
   bowtie2-build ${fasta} ${id}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process bowtie2_align {

   tag "${id}-${genome_acc}" 
   label "big_cpu"

   errorStrategy 'retry'
   maxRetries 1

   publishDir( 
      "${params.outputs}/mapped", 
      mode: 'copy',
      saveAs: { "${id}.${it}" }
   )

   input:
   tuple val( id ), path( reads ), path( idx ), val( genome_acc )

   output:
   tuple val( id ), path( "mapped.sam" ), emit: main
   path "*.log", emit: logs

   script:
   """
   bowtie2 \
      -x ${genome_acc} \
      --very-sensitive-local \
      --threads ${task.cpus} \
      -1 "${reads[0]}" \
      -2 "${reads[1]}" \
      -S mapped.sam \
   2> bowtie2.log

   """
}

