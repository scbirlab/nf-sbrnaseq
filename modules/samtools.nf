process SAMtools_stats {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "stats.txt" )

    script:
    """
    samtools stats -@ ${task.cpus} "${bamfile[0]}" > stats.txt
    """

}

process plot_bamstats {

    tag "${id}"

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( txt )

    output:
    tuple val( id ), path( "stats-*.png" )

    script:
    """
    plot-bamstats -p stats "${txt}"
    """

}



process SAMtools_flagstat {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "flagstat.tsv" )

    script:
    """
    #samtools collate -@ ${task.cpus} -O -u "${bamfile[0]}" \
    #| samtools fixmate -@ ${task.cpus} -m -u - - \
    #| samtools sort -@ ${task.cpus} -m ${Math.round(task.memory.getGiga() * 0.8)}G -u - \
    #| samtools markdup -@ ${task.cpus} - markdup.bam
    #samtools flagstat -@ ${task.cpus} -O tsv markdup.bam > flagstat.tsv
    samtools flagstat -@ ${task.cpus} -O tsv "${bamfile[0]}" > flagstat.tsv
    """

}

process SAMtools_coverage {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "coverage.tsv" )

    script:
    """
    samtools coverage "${bamfile[0]}" -o coverage.tsv
    """

}

process remove_multimappers {

    tag "${id}"
    label 'med_mem'

    // publishDir( 
    //     "${params.outputs}/samtools", 
    //     mode: 'copy',
    //     saveAs: { "${id}.${it}" },
    // )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "sorted.{bam,bai}" )

    script:
    """
    # filter out multimappers
    samtools view -@ ${task.cpus} -F260 -bS -q 3 "${bamfile}" -o filtered.bam
    samtools sort -@ ${task.cpus} -m ${Math.round(task.memory.getGiga() * 0.8)}G filtered.bam -o sorted.bam
    samtools index sorted.bam
    rm filtered.bam
    """

}


