/*
 * Use featureCounts (from the SubReads package) to count transcripts mapping to each gene
 * ## $(ANN_TYPE): which of "CDS", "gene", "mRNA", etc
 * ## $(LABEL): the tag from column 9 to use to label transcript counts, e.g. "Locus", "Name"
 */
process featurecounts {

   tag "${id}"
   label 'big_cpu'

   errorStrategy 'ignore'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/featurecounts", 
      mode: 'copy',
      // pattern: "*.log",
   )

   input:
   tuple val( id ), path( bamfile ), path( gff )
   val nanopore
   val strand
   val annotation_type
   val attribute_label
   val reverse_mate

   output:
   tuple val( id ), path( "*.bam" ), emit: main
   tuple val( id ), path( "*.featureCounts.tsv" ), emit: table
   tuple val( id ), path( "*.species-chr-map.tsv" ), emit: mapping
   path "*.{summary,log}", emit: logs


   script:
   """
   # make species -> chromosome mapping
   awk -F'\\t' -v OFS='\\t' '
      BEGIN { print "taxon_url", "assembly", "genome_accession", "Chr" }
      /^#!genome-build / { split(\$1, _assembly, " "); assembly = _assembly[2]; next } 
      /^#!genome-build-accession / { split(\$1, _acc0, " "); split(_acc0[2], _acc, ":"); acc = _acc[2]; next }
      /^##species / { split(\$1, _url, " "); taxon_url = _url[2]; next }
      !/#/ { chr = \$1; if (!(chr in taxon)) {
            taxon[chr] = taxon_url;
            assem[chr] = assembly;
            genom[chr] = acc;
         }
         next
      } 
      END { for (chr in taxon) print taxon[chr], assem[chr], genom[chr], chr }
      ' \
      "${gff}" \
   > "${id}".species-chr-map.tsv

   samtools view -h ${reverse_mate ? "-f128" : ""} \
      "${bamfile}" -o "${id}.bam"
   samtools index "${id}.bam"
   featureCounts \
      -s ${strand ? strand : "0"} \
      -t ${annotation_type} \
      -g ${attribute_label} \
      -a "${gff}" \
      ${nanopore ? "-L" : "-p"} ${(reverse_mate || nanopore) ? "" : "-B -C --countReadPairs"} -R BAM \
      -T ${task.cpus} \
      -M \
      --verbose \
      --extraAttributes ID,Name,gene_biotype,locus_tag \
      -o "${id}.featureCounts0.tsv" \
      "${id}.bam"
   mv "${id}.featureCounts0.tsv.summary" "${id}.featureCounts.tsv.summary"
   cp "${id}.featureCounts.tsv.summary" "${id}.featureCounts.log"
   rm "${id}.bam"

   python -c '
   import pandas as pd

   (
      pd.read_csv(
         "${id}.featureCounts0.tsv",
         sep="\\t",
         comment="#",
      )
      .merge(
         pd.read_csv("${id}.species-chr-map.tsv", sep="\\t")
      )
      .to_csv("${id}.featureCounts.tsv", sep="\\t", index=False)
   )
   '
   
   """
}
