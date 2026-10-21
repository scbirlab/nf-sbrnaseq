process PlotRNAUMAPs {
    tag "${id}"
    label "some_mem"

    publishDir(
        "${params.outputs}/visualization/rna",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( clustered ), path( h5mu )

    output:
    tuple val( id ), path( "umap.*.png" ), emit: figures
    tuple val( id ), path( "umap.*.csv" ), emit: tables

    script:
    """
    set -euox pipefail

    export MPLCONFIGDIR="mpl"
    export NUMBA_CACHE_DIR="numba"
    mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

    python ${projectDir}/bin/viz/plot-umaps.py \\
        --clustered "${clustered}" \\
        --mudata "${h5mu}" \\
        --output umap

    """
}

process PlotModSites {
    tag "${id}"
    label "some_mem"

    publishDir(
        "${params.outputs}/visualization/modifications",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( clustered ), path( h5mu )

    output:
    tuple val( id ), path( "modsites.*.png" ), emit: figures
    tuple val( id ), path( "modsites.*.csv" ), emit: tables

    script:
    """
    set -euox pipefail

    export MPLBACKEND=Agg
    export MPLCONFIGDIR="mpl"
    export NUMBA_CACHE_DIR="numba"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"
    mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

    python ${projectDir}/bin/viz/plot-mod-sites.py \\
        --clustered "${clustered}" \\
        --mudata "${h5mu}" \\
        --output modsites \\
        --min-depth 1
    """
}


process EmbedTRNA {
    tag "${id}"
    label "med_mem"

    publishDir(
        "${params.outputs}/visualisation/trna",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( clustered )

    output:
    tuple val( id ), path( "trna.h5ad" ), emit: anndata
    tuple val( id ), path( "trna.*.png" ), emit: figures
    tuple val( id ), path( "trna.*.csv" ), emit: tables

    script:
    """
    set -euox pipefail
    export MPLBACKEND=Agg
    export MPLCONFIGDIR=mpl
    export NUMBA_CACHE_DIR=numba
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"
    mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

    python ${projectDir}/bin/viz/embed-trna.py \\
        --clustered "${clustered}" \\
        --output trna
    """
}


process PlotEnd5Usage {
    tag "${id}"
    label "med_mem"

    publishDir("${params.outputs}/visualisation/end5",
        mode: 'copy',
        saveAs: { "${id}.${it}" })

    input:
    tuple val( id ), path( clustered ), path( h5mu )

    output:
    tuple val( id ), path( "end5.*.png" ), emit: figures
    tuple val( id ), path( "end5.*.csv" ), emit: tables

    script:
    """
    set -euox pipefail
    export MPLBACKEND=Agg
    export MPLCONFIGDIR=mpl
    export NUMBA_CACHE_DIR=numba
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"
    mkdir "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

    python ${projectDir}/bin/viz/plot-end5-usage.py \\
        --clustered "${clustered}" \\
        --mudata "${h5mu}" \\
        --output end5
    """
}


process PlotGeneCorrelations {
    tag "${id}"
    label "med_mem"

    publishDir(
        "${params.outputs}/visualization/correlations",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( clustered )

    output:
    tuple val( id ), path( "correlations.png" ), emit: figures

    script:
    """
    set -euox pipefail
    export MPLBACKEND=Agg
    export MPLCONFIGDIR=mpl
    export NUMBA_CACHE_DIR=numba
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p mpl numba "\$XDG_CACHE_HOME/fontconfig"

    python ${projectDir}/bin/viz/plot-gene-correlations.py \\
        --clustered "${clustered}" \\
        --output correlations
    """
}
