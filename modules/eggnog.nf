process DownloadEggnogDBs {

    tag "${url}"

    input:
    val url

    output:
    path "eggnogdb"

    script:
    def _url = url ? url : "http://eggnog6.embl.de/download/emapperdb-5.0.2/"
    """
    set -eoux pipefail
    
    mkdir eggnogdb

    # URL of the directory
    URL="${_url}"

    # Download the directory listing
    curl -v "\$URL" > directory-listing.html

    # Extract file links and download them
    for f in \$(grep -oP '(?<=href=")[^"]*' directory-listing.html | grep -e '^eggnog')
    do
        >&2 echo "[INFO] Downloading \$f..."
        curl -v "\$URL/\$f" -o eggnogdb/\$f
    done

    for f in eggnogdb/*.gz
    do
        gunzip -v "\$f"
    done

    for f in eggnogdb/*.tar
    do
        tar -xvf "\$f"
    done

    for f in eggnog.*
    do
        mv "\$f" eggnogdb/
    done

    """

}

process EggnogFunctionalSets {

    tag "${id}"
    label 'big_cpu_mem'

    publishDir(
        "${params.outputs}/genome", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( fasta ), path( gff ), path( dbs )

    output:
    tuple val( id ), path( "emap.gff" ), emit: gff
    tuple val( id ), path( "emap.tsv" ), emit: table

    script:
    """
    mkdir tmp scratch

    emapper.py \\
        -i "${fasta}" \\
        --data_dir "${dbs}" \\
        --itype proteins \\
        --tax_scope bacteria \\
        --cpu ${task.cpus} \\
        --dbmem \\
        --scratch_dir scratch \\
        --temp_dir tmp \\
        -o emap \\
        --decorate_gff "${gff}"

    mv emap.emapper.annotations emap.tsv
    mv emap.emapper.decorated.gff emap.gff

    """
}
