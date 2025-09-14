// Do quality control checks
process fastQC {

   tag "${sample_id}"
   label 'big_cpu'
   
   input:
   tuple val( sample_id ), path ( reads )

   output:
   tuple val( sample_id ), path ( "*.zip" ), emit: logs
   path "*.zip", emit: multiqc_logs

   script:
   """
   zcat ${reads} > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """
   stub:
   """
   zcat ${reads} | head -n1000 > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """

}