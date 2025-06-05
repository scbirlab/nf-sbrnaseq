/*
 * Make log report
 */
process multiQC {

   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/multiqc", 
      mode: 'copy',
   )

   container 'docker://multiqc/multiqc:v1.29'

   input:
   path( '*', stageAs: '?/*' )

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}