
process CountModSites {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/modifications",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bam ), path( bai ), path( sites ), path( raw_bam )
    val min_mapq
    val min_baseq
    val mapper
    val mapping_min_mapq
    val allow_multimap

    output:
    tuple val( id ), path( "mod-sites.tsv.gz" )

    script:
    """
    set -euox pipefail

    samtools sort \\
        -@ ${task.cpus} \\
        -m ${Math.round(Math.floor(task.memory.getMega() * 0.8 / task.cpus))}M \\
        -o raw.sorted.bam \\
        "${raw_bam}"

    samtools index raw.sorted.bam

    python ${projectDir}/bin/count-mod-sites.py \\
        --bam "${bam}" \\
        --raw-bam raw.sorted.bam \\
        --sites "${sites}" \\
        --sample-id "${id}" \\
        --min-mapq ${min_mapq} \\
        --min-baseq ${min_baseq} \\
        --mapper "${mapper}" \\
        --mapping-min-mapq ${mapping_min_mapq} \\
        --allow-multimap ${allow_multimap} \\
        --output "mod-sites.tsv.gz"

    rm raw.sorted.bam{,.bai}
    """
}
