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

    container 'https://depot.galaxyproject.org/singularity/rseqc%3A5.0.4--pyhdfd78af_0'

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

    container 'https://depot.galaxyproject.org/singularity/bedops%3A2.4.42--h9948957_0'

    script:
    """
    gff2bed < ${gff} > bedfile.bed
    """
}