process build_AnnData {
   tag "${id}"
   label "big_mem"

   publishDir( 
      "${params.outputs}/scanpy", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( counts_table )

   output:
   tuple val( id ), path( "*.h5ad" ), emit: main
   tuple val( id ), path( "*.{png,csv}" ), emit: plots
   path "*.log", emit: logs

   script:
   """
   export MPLCONFIGDIR="mpl"
   export NUMBA_CACHE_DIR="numba"
   mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

   python ${projectDir}/bin/anndata-utils.py build "${counts_table}" --id "${id}" -o "build" 2> build.log

   """
}

process filter_AnnData {
   tag "${id}"
   label "med_mem"

   publishDir( 
      "${params.outputs}/scanpy", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( anndata )

   output:
   tuple val( id ), path( "*.h5ad" ), emit: main
   tuple val( id ), path( "*.{png,csv}" ), emit: plots, optional: true
   path "*.log", emit: logs

   script:
   """
   export MPLCONFIGDIR="mpl"
   export NUMBA_CACHE_DIR="numba"
   mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

   python ${projectDir}/bin/anndata-utils.py filter "${anndata}" -o filter 2> filter.log

   """
}


process cluster_cells {

   tag "${id}"
   label "med_mem"

   publishDir( 
      "${params.outputs}/clustering", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( anndata )

   output:
   tuple val( id ), path( "*.{png,csv}" ), emit: main
   path "*.log", emit: logs

   script:
   """
   export MPLCONFIGDIR="mpl"
   export NUMBA_CACHE_DIR="numba"
   mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

   python ${projectDir}/bin/anndata-utils.py cluster "${anndata}" -o cluster 2> cluster.log
   mv figures/*.png .
   
   """
}