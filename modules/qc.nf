// Do quality control checks
process fastQC {

   tag "${id}"
   label 'big_cpu'
   
   input:
   tuple val( id ), path ( reads, stageAs: '??/*' )

   output:
   tuple val( id ), path ( "*.zip" ), emit: logs
   path "*.zip", emit: multiqc_logs

   script:
   """
   zcat ??/*.fastq.gz > "${id}.fastq"
   fastqc --noextract --memory 10000 --threads ${task.cpus} "${id}.fastq"
   rm "${id}.fastq"
   """
   stub:
   """
   zcat ??/*.fastq.gz | head -n1000 > "${id}.fastq"
   fastqc --noextract --memory 10000 --threads ${task.cpus} "${id}.fastq"
   rm "${id}.fastq"
   """

}