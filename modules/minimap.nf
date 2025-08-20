process minimap_index {

   tag "${genome_id}"
   
   input:
   tuple val( genome_id ), path( fasta )
   val nanopore

   output:
   tuple val( genome_id ), path( "target.mmi" )

   script:
   """
   minimap2 -x ${nanopore ? 'map-ont' : 'sr'} -d target.mmi ${fasta}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process minimap_align {

   tag "${id}:${genome_acc}" 
   label "big_mem"
   time "2d"

   // errorStrategy 'retry'
   // maxRetries 1

   publishDir( 
      "${params.outputs}/mapped", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
      pattern: "*.{sam,log}"
   )

   input:
   tuple val( id ), path( reads ), path( idx ), val( genome_acc )
   val nanopore

   output:
   tuple val( id ), path( "minimap2.sam" ), emit: main
   path "minimap2.log", emit: logs

   script:
   """
   minimap2 -a -c -2 --MD -x ${nanopore ? 'map-ont' : 'sr --frag=yes'} \
      ${idx} ${reads} \
      -o minimap2.sam \
   2> minimap2.log

   """
}


