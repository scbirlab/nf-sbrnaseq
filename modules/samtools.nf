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

    input:
    tuple val( id ), path( bamfile )
    val mapper
    val min_mapq
    val require_unique

    output:
    tuple val( id ), path( "*.sorted.{bam,bai}" )

    script:
    def unique = (require_unique && mapper == 'star' ? "-e '[NH]==1'" : "")
    """
    set -euox pipefail

    # Exclude unmapped, secondary and supplementary alignments.
    # Apply minimum mapping confidence.
    #
    # STAR: NH==1 is required in strict mode.
    # Bowtie2/minimap2: use MAPQ and primary status;
    # neither provides a universally reliable NH tag.

    samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -F ${require_unique ? '2308' : '260'} \\
        -q ${min_mapq} ${unique} \\
        "${bamfile}" \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -m ${Math.round(Math.floor(task.memory.getMega() * 0.8 / task.cpus))}M \\
        -o "${id}.sorted.bam" -

    samtools index "${id}.sorted.bam" -o "${id}.sorted.bai"

    samtools quickcheck -v "${id}.sorted.bam"

    """
}
