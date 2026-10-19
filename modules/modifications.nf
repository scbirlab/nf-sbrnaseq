
process CountModSites {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/modifications",
        mode: 'copy',
    )

    input:
    tuple val( id ), path( bam ), path( bai ), path( sites )
    val min_mapq
    val min_baseq

    output:
    tuple val( id ), path( "${id}.mod-sites.tsv.gz" )

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/count-mod-sites.py \\
        --bam "${bam}" \\
        --sites "${sites}" \\
        --sample-id "${id}" \\
        --min-mapq ${min_mapq} \\
        --min-baseq ${min_baseq} \\
        --output "${id}.mod-sites.tsv.gz"
    """
}
