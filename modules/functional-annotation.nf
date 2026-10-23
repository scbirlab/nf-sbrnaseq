
process EmapperGeneTable {
    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/genome",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val(id), path(gff), path(emapper)

    output:
    tuple val(id), path("functional-annotations.tsv"), emit: table

    script:
    """
    set -euox pipefail

    python ${projectDir}/bin/emapper-gene-table.py \\
        --gff "${gff}" \\
        --emapper "${emapper}" \\
        --output functional-annotations.tsv
    """
}
