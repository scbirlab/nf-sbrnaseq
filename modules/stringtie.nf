process stringtie {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/stringtie", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile ), path( gff )

    output:
    tuple val( id ), path( "stringtie.gtf" ), emit: transcripts
    tuple val( id ), path( "{gene_abund.tab,cov_refs.gtf}" ), emit: coverage

    script:
    """
    stringtie \
        -o stringtie.gtf \
        -A gene_abund.tab \
        -C cov_refs.gtf \
        -p ${task.cpus} \
        --fr -m 50 \
        -G "${gff}" \
        "${bamfile[0]}" 
    """

}

process merge_stringtie {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/stringtie", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( gtfs ), path( gff )

    output:
    tuple val( id ), path( "st-merge.gtf" )

    script:
    """
    stringtie --merge -o st-merge.gtf -G "${gff}" ${gtfs}
    """

}


process stringtie_count {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/stringtie", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile ), path( gff )

    output:
    tuple val( id ), path( "stringtie-count.gtf" ), emit: gtf
    tuple val( id ), path( "{gene_abund-count.tab,cov_refs-count.gtf}" ), emit: coverage
    tuple val( id ), path( "*.ctab" ), emit: counts


    script:
    """
    stringtie \
        -e -B \
        -o stringtie-count.gtf \
        -A gene_abund-count.tab \
        -C cov_refs-count.gtf \
        -p ${task.cpus} \
        --fr -m 50 \
        -G "${gff}" \
        "${bamfile[0]}" 
    """

}