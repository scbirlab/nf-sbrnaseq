
process ExtractMoleculeEnds {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/ends",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bam ), path( bai )
    val read_number
    val boundary
    val strandedness
    val min_mapq
    val exclude_softclip

    output:
    tuple val( id ), path( "ends.tsv.gz" )

    script:
    """
    set -euox pipefail

    python ${projectDir}/bin/extract-ends.py \\
        --bam "${bam}" \\
        --sample-id "${id}" \\
        --read ${read_number} \\
        --boundary "${boundary}" \\
        --strand ${strandedness} \\
        --min-mapq ${min_mapq} \
        ${exclude_softclip ? '--exclude-softclip' : ''} \\
        --output "ends.tsv.gz"
    
    """
}


process ClusterMoleculeEnds {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/ends",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( ends, stageAs: 'ends/?????.tsv' )
    val min_support
    val max_width

    output:
    tuple val( id ), path( "loci.tsv" ), path( "loci.bed" ), emit: loci
    tuple val( id ), path( "paraclu.tsv" ), emit: raw

    script:
    """
    set -euox pipefail

    zcat ${ends} \\
    | awk -F '\\t' -v OFS='\\t' '
        NR == 1 {
            for ( i=1; i<=NF; i++ ) header[\$i] = i
        }
        NR > 1 {
            chrom = \$header["contig"]
            strand = \$header["strand"]
            pos = \$header["pos1"]
            counter[chrom, strand, pos]++
        }
        END {
            for (key in counter) {
                split(key, a, SUBSEP)
                print a[1], a[2], a[3], counter[key]
            }
        }
    ' \\
    | grep -v '^contig' \\
    | sort -k1,1 -k2,2 -k3,3n \\
    > paraclu-input.tsv

    if [ "\$(wc -l < paraclu-input.tsv)" -gt 0 ]
    then
        paraclu "${min_support}" paraclu-input.tsv > "paraclu.tsv"
    else
        touch "paraclu.tsv"
    fi

    python ${projectDir}/bin/filter-ends.py \\
        "paraclu.tsv" \\
        --max-width ${max_width}

    sort -k1,1 -k2,2n "loci.bed" > "loci.sorted.bed"
    mv "loci.sorted.bed" "loci.bed"
    
    """
}


process AnnotateEndLoci {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/ends",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( loci_tsv ), path( loci_bed ), path( gff )

    output:
    tuple val( id ), path( "annotated.tsv" )

    script:
    """
    set -euox pipefail

    grep -v '^#' "${gff}" \\
    | awk -v OFS='\\t' \\
        '\$3 == "gene" && NF == 9 &&
        (\$7 == "+" || \$7 == "-") {
            print \$1, \$4-1, \$5, \$9, \$6, \$7
        }' \\
    > "features.bed"

    bedtools closest \\
        -a "${loci_bed}" \\
        -b "features.bed" \\
        -s -D a -t first \\
    > annotated.bed

    python ${projectDir}/bin/annotate-end-loci.py \\
        annotated.bed \\
        --output annotated.tsv
    
    """
}



process AssignEndLoci {

    tag "${id}"
    label 'some_mem'

    publishDir(
        "${params.outputs}/ends",
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( ends, stageAs: 'ends-??????/*' ), path( loci )

    output:
    tuple val( id ), path( "counts.tsv" ), emit: counts
    tuple val( id ), path( "assignment-qc.tsv" ), emit: qc

    script:
    """
    set -euox pipefail

    python ${projectDir}/bin/assign-end-loci.py \\
        --loci "${loci}" \\
        --ends ends-*/*.gz \\
        --counts counts.tsv \\
        --qc assignment-qc.tsv
    """
}

