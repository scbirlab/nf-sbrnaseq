
process ExtractMoleculeEnds {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/ends",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bam ), path( bai )
    val read_number
    val boundary
    val strandedness
    val min_mapq
    val exclude_softclip

    output:
    tuple val( id ), path( "ends.tsv.gz" )

    script:
    """
    set -euox pipefail

    python ${projectDir}/bin/extract-ends.py \\
        --bam "${bam}" \\
        --sample-id "${id}" \\
        --read ${read_number} \\
        --boundary "${boundary}" \\
        --strand ${strandedness} \\
        --min-mapq ${min_mapq} \
        ${exclude_softclip ? '--exclude-softclip' : ''} \\
        --output "ends.tsv.gz"
    
    """
}
