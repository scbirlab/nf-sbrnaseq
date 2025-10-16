process SAMtools_stats {

    tag "${id}"
    label 'big_cpu'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        // saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.stats.txt" )

    script:
    """
    samtools stats -@ ${task.cpus} "${bamfile[0]}" > "${id}.stats.txt"
    """

}

process plot_bamstats {

    tag "${id}"

    errorStrategy 'ignore'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        // saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( txt )

    output:
    tuple val( id ), path( "${id}-*.png" )

    script:
    """
    plot-bamstats -p "${id}" "${txt}"
    """

}


process SAMtools_flagstat {

    tag "${id}"
    label 'big_cpu'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        // saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.flagstat.tsv" )

    script:
    """
    #samtools collate -@ ${task.cpus} -O -u "${bamfile[0]}" \
    #| samtools fixmate -@ ${task.cpus} -m -u - - \
    #| samtools sort -@ ${task.cpus} -m ${Math.round(task.memory.getGiga() * 0.8)}G -u - \
    #| samtools markdup -@ ${task.cpus} - markdup.bam
    #samtools flagstat -@ ${task.cpus} -O tsv markdup.bam > flagstat.tsv
    samtools flagstat -@ ${task.cpus} -O tsv "${bamfile[0]}" > "${id}.flagstat.tsv"
    """

}

process SAMtools_coverage {

    tag "${id}"

    publishDir( 
        "${params.outputs}/samtools-coverage", 
        mode: 'copy',
        // saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.tsv" )

    script:
    """
    samtools coverage "${bamfile[0]}" -o "${id}.tsv"
    """

}

process remove_multimappers {

    tag "${id}"
    label 'big_mem'

    // publishDir( 
    //     "${params.outputs}/samtools", 
    //     mode: 'copy',
    //     saveAs: { "${id}.${it}" },
    // )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "*.sorted.{bam,bai}" )

    script:
    """
    # filter out multimappers
    #samtools view -@ ${task.cpus} -F2308 -bS --min-MQ=1 "${bamfile}" -o filtered.bam
    samtools view -@ ${task.cpus} -bS "${bamfile}" -o filtered.bam
    samtools sort -@ ${task.cpus} -m 2G filtered.bam -o "${id}.sorted.bam"
    samtools index "${id}.sorted.bam" -o "${id}.sorted.bai"
    rm filtered.bam
    
    """

}
