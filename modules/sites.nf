
process ResolveSites {

    tag "${id}"

    publishDir(
        "${params.outputs}/modifications",
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( fasta ), path( gff )
    path rules
    val join_on

    output:
    tuple val( id ), path( "resolved-sites.tsv" )

    script:
    """
    set -euox pipefail

    samtools faidx "${fasta}"

    python ${projectDir}/bin/resolve-sites.py \\
        --gff "${gff}" \\
        --fasta "${fasta}" \\
        --sites "${rules}" \\
        --join-on "${join_on}" \\
        --output resolved-sites.tsv
    
    """
}
