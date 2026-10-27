
process BuildEnd5AnnData {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/mudata",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( counts ), path( qc ), path( loci )

    output:
    tuple val( id ), path( "end5.h5ad" )

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/anndata/build-end5-anndata.py \\
        --counts "${counts}" \\
        --qc "${qc}" \\
        --loci "${loci}" \\
        --sample-id "${id}" \\
        --output end5.h5ad
    """
}


process BuildModAnnData {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/mudata",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( counts )

    output:
    tuple val( id ), path( "mod.h5ad" )

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/anndata/build-mod-anndata.py \\
        --counts "${counts}" \\
        --sample-id "${id}" \\
        --output mod.h5ad
    """
}


process BuildMuData {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/mudata",
        mode: 'copy',
        saveAs: { "${id}.h5mu" },
    )

    input:
    tuple val( id ), path( rna ), path( rna_qc ), path( h5ads, stageAs: "modalities/*" )

    output:
    tuple val( id ), path( "sample.h5mu" )

    script:
    """
    set -euox pipefail

    python ${projectDir}/bin/anndata/build-mudata.py \\
        --rna "${rna}" \\
        --rna-qc "${rna_qc}" \\
        --modalities modalities/*.h5ad \\
        --sample-id "${id}" \\
        --output sample.h5mu
    """
}


process QCMuData {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/mudata/qc",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( h5mu ) 

    output:
    tuple val( id ), path( "qc.*.png" ), emit: figures
    tuple val( id ), path( "qc.*.csv" ), emit: tables

    script:
    """
    set -euox pipefail

    export MPLBACKEND=Agg
    export MPLCONFIGDIR=mpl
    export NUMBA_CACHE_DIR=numba
    mkdir -p mpl numba

    python ${projectDir}/bin/anndata/qc-mudata.py \\
        "${h5mu}" \\
        --output qc
    """
}
