process gene_body_coverage {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/coverage", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile ), path( bed )

    output:
    tuple val( id ), path( "*.pdf" ), emit: main
    path "*.geneBodyCoverage.txt", emit: logs

    script:
    """
    geneBody_coverage.py -r "${bed}" -i "${bamfile}" -o gene-body
    """

}

process gff2bed {
    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/coverage", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( gff )

    output:
    tuple val( id ), path( "*.bed" )

    script:
    """
    gff2bed < ${gff} > bedfile.bed
    """
}