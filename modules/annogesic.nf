process annogesic_transcript {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/transcripts", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    container 'silasysh/annogesic:latest'

    input:
    tuple val( id ), path( wig )

    output:
    tuple val( id ), path( "*.gff" ), emit: gff
    tuple val( id ), path( "{stat_,}*.{csv,png}" ), emit: stats
    tuple val( id ), path( "*_transcript.csv" ), emit: table

    script:
    """
    annogesic create --project_path ANNOgesic
    annogesic transcript --project_path ANNOgesic --frag_libs "${wig}"
    
    mv ANNOgesic/output/transcripts/{gffs/*.gff,tables/*.csv,statistics/*.png,statistics/*.csv} .
    
    """
}

process bam2wig {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/transcripts", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bam )

    output:
    tuple val( id ), path( "coverage.bg" )

    script:
    """
    bamCoverage \
        -b "${bam}" \
        --outFileName coverage.bg \
        --outFileFormat bedgraph \
        --binSize 1 \
        --numberOfProcessors ${task.cpus}
    
    """
}
