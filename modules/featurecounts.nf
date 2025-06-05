/*
 * Use featureCounts (from the SubReads package) to count transcripts mapping to each gene
 * ## $(ANN_TYPE): which of "CDS", "gene", "mRNA", etc
 * ## $(LABEL): the tag from column 9 to use to label transcript counts, e.g. "Locus", "Name"
 */
process featurecounts {

   tag "${id}"
   label 'med_mem'

   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      pattern: "*.log",
   )

   container 'docker://genomicpariscentre/featurecounts:1.5.3'

   input:
   tuple val( id ), path( bamfile ), path( gff )
   val paired

   output:
   tuple val( id ), path( "*.bam" ), emit: main
   tuple val( id ), path( "*.featureCounts.tsv" ), emit: table
   path "*.summary", emit: logs

   script:
   """
   # make species -> chromosome mapping
   cat \
      <(printf 'genome_accession\\tChr\\n') \
      <(awk -v OFS='\\t' \
         '
         /^#!genome-build-accession/ { split(\$2, _species, ":"); species = _species[2] } 
         !/#/ { if (!(\$1 in chr)) chr[\$1] = species } 
         END { for (key in chr) print chr[key], key}
         ' \
         "${gff}" \
      | sort -k2) \
   > species-chr-map.tsv

   samtools index "${bamfile}"
   featureCounts \
      -p -C \
      -s ${params.strand} \
      -t ${params.ann_type} \
      -g ${params.label} \
      -a "${gff}" \
      -R BAM \
      -T ${task.cpus} ${paired ? '--countReadPairs' : ''} \
      --verbose \
      --extraAttributes Name,gene_biotype,locus_tag \
      -o "${id}.featureCounts0.tsv" \
      "${bamfile}"
   mv "${id}.featureCounts0.tsv.summary" "${id}.featureCounts.tsv.summary"
   cp "${id}.featureCounts.tsv.summary" "${id}.featureCounts.log"

   cat \
      <(head -n1 "${id}.featureCounts0.tsv") \
      <(cat \
         <(head -n2 "${id}.featureCounts0.tsv" | tail -n1) \
         <(tail -n+3 "${id}.featureCounts0.tsv" | sort -k2) \
       | join species-chr-map.tsv - -j2 --header -t\$'\\t') \
   > "${id}.featureCounts.tsv"
   
   """
}
