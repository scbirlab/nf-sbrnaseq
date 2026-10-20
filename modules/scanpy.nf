process build_AnnData {
   tag "${id}"
   label "big_mem"

   errorStrategy 'ignore'

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
   path "build.log", emit: logs

   script:
   """
   set -euox pipefail

   export MPLCONFIGDIR="mpl"
   export NUMBA_CACHE_DIR="numba"
   mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

   python ${projectDir}/bin/anndata-utils.py build "${counts_table}" --id "${id}" -o "build"  2> >(tee build.log)

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
   set -euox pipefail

   export MPLCONFIGDIR="mpl"
   export NUMBA_CACHE_DIR="numba"
   mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

   python ${projectDir}/bin/anndata-utils.py filter "${anndata}" \\
      --min-counts-per-cell ${params.min_counts_per_cell} \\
      --min-cells-per-gene ${params.min_cells_per_gene} \\
      -o filter \\
   2> >(tee filter.log)


   """
}


process cluster_cells {

   tag "${id}"
   label "med_mem"

   errorStrategy 'ignore'

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
   set -euox pipefail
   
   export MPLCONFIGDIR="mpl"
   export NUMBA_CACHE_DIR="numba"
   mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

   python ${projectDir}/bin/anndata-utils.py cluster "${anndata}" -o cluster 2> >(tee cluster.log)

   mv figures/*.png .
   
   """
}