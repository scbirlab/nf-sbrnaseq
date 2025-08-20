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
    set -x
    samtools index "${bamfile}" 
    geneBody_coverage.py -r "${bed}" -i "${bamfile}" -o gene-body

    output_lines=\$(wc -l < *.geneBodyCoverage.txt)
    if [ "\$output_lines" -eq 1 ]
    then
        echo "Failed!"
        exit 1
    fi
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
    gff2bed < "${gff}" \
    | awk -v OFS='\\t' '\$8 == "gene" { \$7=\$2; \$8=\$3; \$9="0,0,0"; \$10=1; \$11=\$3-\$2","; \$12="0,"; print }' \
    > bedfile.bed
    
    """
}