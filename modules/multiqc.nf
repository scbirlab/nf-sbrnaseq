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

   input:
   path( '*', stageAs: '?/*' )

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}