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
      "${params.outputs}/mapped/bowtie2", 
      mode: 'copy',
      // saveAs: { "${id}.${it}" }
   )

   input:
   tuple val( id ), path( reads ), path( idx ), val( genome_acc )
   val max_multimaps

   output:
   tuple val( id ), path( "*.bam" ), emit: main
   path "*.log", emit: logs

   script:
   """
   bowtie2 \\
      -x ${genome_acc} \\
      --very-sensitive-local \\
      -k ${max_multimaps} \\
      --threads ${task.cpus} \\
      --no-unal \\
      -1 "${reads[0]}" \\
      -2 "${reads[1]}" \\
      -S mapped.sam \\
   2> bowtie2.log

   samtools view -bS "mapped.sam" -o "${id}.bowtie2.bam"
   rm "mapped.sam"

   """
}

