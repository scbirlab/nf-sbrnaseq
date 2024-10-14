#!/usr/bin/env nextflow

/*
========================================================================================
   Single Bacterium RNASeq Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-sbrnaseq
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
pipeline_title = """\
                 S C B I R   s b R N A - S E Q   P I P E L I N E
                 ===============================================
                 """
                 .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println pipeline_title + """\
         Nextflow pipeline to process demultiplexed Illumina paired-end 
         FASTQ files from multiple bacterial samples into a gene x
         cell count table.

         Usage:
            nextflow run sbcirlab/nf-sbrnaseq --sample_sheet <csv> --fastq-dir <dir>
            nextflow run sbcirlab/nf-sbrnaseq -c <config-file>

         Required parameters:
            sample_sheet         Path to a CSV containing sample IDs matched with FASTQ filenames, genome information, and adapter sequences.
            fastq_dir            Path to directory containing the FASTQ file.

         Optional parameters (with defaults):   
            trim_qual = 5        For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = "9:38"  For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
            strand = 1           For `featureCounts`, the strandedness of RNA-seq. `1` for forward, `2` for reverse.
            ann_type = 'gene'    For `featureCounts`, features from GFF column 3 to use for counting
            label = 'Name'       For `featureCounts`, one or more (comma-separated) fields from column 9 of GFF for labeling counts

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   System.exit(0)
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}
if ( !params.fastq_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to fastq_dir")
}

working_dir = params.outputs

processed_o = "${working_dir}/processed"
counts_o = "${working_dir}/counts"
multiqc_o = "${working_dir}/multi_qc"

log.info pipeline_title + """\
   inputs
      input dir.     : ${params.inputs}
      FASTQ dir.     : ${params.fastq_dir}
      sample sheet   : ${params.sample_sheet}
   trimming 
      quality        : ${params.trim_qual}
      minimum length : ${params.min_length}
   output            
      Processed      : ${processed_o}
      Counts         : ${counts_o}
      MultiQC        : ${multiqc_o}
   """
   .stripIndent()

dirs_to_make = [processed_o, counts_o, multiqc_o]

log.info  """\
          Making directories: 
          """
          .stripIndent()

dirs_to_make.each { 
   log.info "${it}: " 
   log.info file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   Channel.of( params.umicollapse_repo )
   | DOWNLOAD_UMICOLLAPSE

   Channel.fromPath( 
         params.sample_sheet, 
         checkIfExists: true 
      )
      .splitCsv( header: true )
      .set { csv_ch }
   csv_ch
      .map { tuple( 
         it.sample_id,
         tuple( it.adapter_read1_5prime, it.adapter_read2_5prime ),
         tuple( it.adapter_read1_3prime, it.adapter_read2_3prime ) 
      ) }
      .set { adapter_ch }  // sample_name, [adapt5], [adapt3] 
   csv_ch
      .map { tuple( 
         it.sample_id,
         tuple( it.umi_read1, it.umi_read2 )
      ) }
      .set { umi_ch }  // sample_id, [umis]
   csv_ch
      .map { tuple( 
         it.sample_id,
         tuple( 
            it.bc1, 
            it.bc2, 
            it.bc3
         ).collect {
            file( 
               "${params.inputs}/${it}",
               checkIfExists: true 
            )
         }
      ) }
      .map { tuple(
         it[0], 
         it[1].collect { x -> x.getSimpleName() }.join("-"),
         it[1]
      ) }
      .set { barcode_ch }  // sample_id, wl_id, [BCs]
   csv_ch
      .map { tuple( 
         it.sample_id,
         file( 
            "${params.fastq_dir}/*${it.fastq_pattern}*",
             checkIfExists: true 
         ).sort()
      ) }
      .set { reads_ch }  // sample_id, [reads]
   csv_ch
      .map { tuple( 
         it.sample_id,
         it.genome_accession 
      ) }
      .set { genome_ch }  // sample_id, genome_acc

   genome_ch
      .map { it[1] }  // genome_acc
      .unique()
      | DOWNLOAD_GENOME   // genome_acc, genome, gff
   DOWNLOAD_GENOME.out
      .map { it[0..1] }  // genome_acc, genome
      | BOWTIE2_INDEX   // genome_acc, [genome_idx]
   genome_ch
      .map { it[1..0] }  // genome_acc, sample_id
      .combine( BOWTIE2_INDEX.out, by: 0 )  // genome_acc, sample_id, [genome_idx]
      .map { it[1..-1] + [ it[0] ] }  // sample_id, [genome_idx], genome_acc
      .set { genome_idx }
   genome_ch
      .map { it[1..0] }  // genome_acc, sample_id
      .combine( DOWNLOAD_GENOME.out, by: 0 )  // genome_acc, sample_id, genome, gff
      .map { tuple( it[1], it[-1] ) }  // sample_id, gff
      .set { genome_gff }
   
   reads_ch | FASTQC 

   Channel.value( 
      tuple( params.trim_qual, params.min_length )
   ).set { trim_params } 
   TRIM_CUTADAPT(
      reads_ch.combine( adapter_ch, by: 0 ),  // sample_id, [reads], [adapt5], [adapt3]
      trim_params
   )  // sample_id, [reads]
   barcode_ch
      .map { it[1..2] }  // wl_id, [BCs]
      .unique()
      | BUILD_WHITELIST  // wl_id, whitelist
      | ADD_WHITELIST_ERRORS // wl_id, whitelist_err
   barcode_ch
      .map { it[1..0] }  // wl_id, sample_id
      .combine( ADD_WHITELIST_ERRORS.out, by: 0 )  // wl_id, sample_id, whitelist_err
      .map { it[1..-1] }  // sample_id, whitelist_err
      .set { whitelists }
   // TRIM_CUTADAPT.out.main | UMITOOLS_WHITELIST
   TRIM_CUTADAPT.out.main  // sample_id, [reads]
      .combine( umi_ch, by: 0 )  // sample_id, [reads], umis
      .combine( whitelists, by: 0 )  // sample_id, [reads], umis, whitelist
      | UMITOOLS_EXTRACT  // sample_id, [reads]

   UMITOOLS_EXTRACT.out.main
      .combine( genome_idx, by: 0 )  // sample_id, [reads], [genome_idx], genome_acc
      | BOWTIE2_ALIGN  // sample_id, [bam_bai]

   // BOWTIE2_ALIGN.out.main
   //    | UMITOOLS_DEDUP  // sample_id, dedup_bam

   BOWTIE2_ALIGN.out.main
      .combine( DOWNLOAD_UMICOLLAPSE.out )  // sample_id, [bam_bai], umicollapse_repo
      | UMICOLLAPSE  // sample_id, dedup_bam

   UMICOLLAPSE.out.main
      .combine( genome_gff, by: 0 )  // sample_id, dedup_bam, gff
      | FEATURECOUNTS  // sample_id, [count_bam_bai]

   FEATURECOUNTS.out.main
      | UMITOOLS_COUNT   // sample_id, counts

   TRIM_CUTADAPT.out.logs
      .concat(
         FASTQC.out.multiqc_logs,
         BOWTIE2_ALIGN.out.logs,
         FEATURECOUNTS.out.logs 
      )
      .flatten()
      .unique()
      .collect()
      | MULTIQC

}

process DOWNLOAD_UMICOLLAPSE {

   input:
   val umicollapse_repo

   output:
   path "UMICollapse"

   script:
   """
   git clone ${umicollapse_repo} #https://github.com/siddharthab/UMICollapse.git
   UMI_LIB_PATH=UMICollapse/lib
   mkdir -p \$UMI_LIB_PATH
   curl -L https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar > \$UMI_LIB_PATH/htsjdk-2.19.0.jar
   curl -L https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar > \$UMI_LIB_PATH/snappy-java-1.1.7.3.jar
   """


}

process DOWNLOAD_GENOME {

   tag "${accession}"
   label 'some_mem'

   input:
   val accession

   output:
   tuple val( accession ), path( "ncbi_dataset/data/*/${accession}_*_genomic.fna" ), path( "ncbi_dataset/data/*/*.gff" )

   script:
   """
   wget "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${accession}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" -O genome-out
   unzip -o genome-out ncbi_dataset/data/${accession}/{${accession}_*_genomic.fna,*.gff}
   """
}


// Do quality control checks
process FASTQC {

   label 'med_mem'

   tag "${sample_id}"

   input:
   tuple val( sample_id ), path ( reads )

   output:
   tuple val( sample_id ), path ( "*.zip" ), emit: logs
   path "*.zip", emit: multiqc_logs

   script:
   """
   zcat ${reads} > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """
   stub:
   """
   zcat ${reads} | head -n1000 > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """

}

/*
 * Trim adapters from reads
 */
process TRIM_CUTADAPT {

   tag "${sample_id}"
   
   errorStrategy 'retry'
   maxRetries 2

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads, stageAs: "?/*" ), val( adapters5 ), val( adapters3 )
   tuple val( trim_qual ), val( min_length )

   output:
   tuple val( sample_id ), path( "*.with-adapters.fastq.gz" ), emit: main
   path( "*.log" ), emit: logs

   script:
   def adapter_5prime_R = (adapters5[0].length() > 0 ? "-G '${adapters5[0]}'" : "")
   """
   for i in \$(seq 1 2)
   do
      zcat */*_R"\$i"*.fastq.gz | gzip --best > ${sample_id}_polyA_R"\$i".fastq.gz
   done

   cutadapt \
		-a "A{10}" \
      -A "T{10}" \
		-q ${trim_qual} \
		--nextseq-trim ${trim_qual} \
		--minimum-length ${min_length} \
		--report full \
      --action trim \
		-o ${sample_id}_3prime_R1.fastq.gz \
      -p ${sample_id}_3prime_R2.fastq.gz \
		${sample_id}_polyA_R?.fastq.gz > ${sample_id}.polyA.cutadapt.log

   cutadapt \
		-a "${adapters3[0]}" \
      -A "${adapters3[1]}" \
      --minimum-length ${min_length} \
		--report full \
      --action trim \
		-o ${sample_id}_5prime_R1.fastq.gz \
      -p ${sample_id}_5prime_R2.fastq.gz \
		${sample_id}_3prime_R?.fastq.gz > ${sample_id}.3prime.cutadapt.log

   cutadapt \
		-g "${adapters5[1]}" ${adapter_5prime_R} \
		--report full \
      --action retain \
		-o ${sample_id}_R1.with-adapters.fastq.gz \
      -p ${sample_id}_R2.with-adapters.fastq.gz \
		${sample_id}_5prime_R?.fastq.gz > ${sample_id}.5prime.cutadapt.log

   rm ${sample_id}_polyA_R?.fastq.gz
   """
}

/*
 * Identify cell barcodes
 */
process UMITOOLS_WHITELIST {
   errorStrategy 'ignore'

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( genome_id ), val( sample_id ), val( adapters ), val( umis ), path( reads )

   output:
   tuple path( "*.txt" ), val( sample_id ), emit: main
   path( "*.log" ), emit: logs
   path( "*.png" ), emit: plots

   script:
   """
   umi_tools whitelist \
      --knee-method distance \
      --allow-threshold-error \
      --plot-prefix ${sample_id}.whitelist \
      --method umis \
		--bc-pattern "${umis[0]}" \
		--bc-pattern2 "${umis[1]}" \
      --extract-method regex \
      --error-correct-threshold ${params.umitools_error} \
      --ed-above-threshold correct \
		--log ${sample_id}.whitelist.log \
		--stdin ${reads[0]} \
		--read2-in ${reads[1]} \
      --stdout ${sample_id}.whitelist.txt
   """
}

// Build whitelist from provided barcode files
process BUILD_WHITELIST {
   tag "${whitelist_id}"

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.txt" )

   input:
   tuple val( whitelist_id ), path( bcs )

   output:
   tuple val( whitelist_id ), path( "*.whitelist.txt" )

   script:
   """
   for f in ${bcs}
   do
      tail -n+2 \$f \
      | tr -d \$'\\r' \
      | sort \
      | awk -F, -v OFS=, '{ print "__JOIN__", \$2 }' \
      | sort -t, -k1 \
      > \$f.temp
   done

   join -t, ${bcs[0]}.temp ${bcs[1]}.temp \
   | join -t, - ${bcs[2]}.temp \
   | awk -F, '{ print \$2\$3\$4 }' \
   > ${whitelist_id}.whitelist.txt
   """
}

// Build whitelist from provided barcode files
process ADD_WHITELIST_ERRORS {
   tag "${whitelist_id}"

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.txt" )

   input:
   tuple val( whitelist_id ), path( whitelist )

   output:
   tuple val( whitelist_id ), path( "*.whitelist-err.txt" )

   script:
   """
   #!/usr/bin/env python

   from itertools import product

   ALPHABET = "ATCG"
   INPUT_FILE = "${whitelist}"
   OUTPUT_FILE = "${whitelist_id}.whitelist-err.txt"

   with open(OUTPUT_FILE, 'w') as outfile, open(INPUT_FILE, 'r') as infile:
      for line in infile:
         line = line.strip()
         alternatives = (
            f"{line[:i]}{letter}{line[(i+1):]}" 
            for (i, char), letter in product(
               enumerate(line), 
               ALPHABET,
            ) if char != letter
         )
         alternatives = ','.join(sorted(set(alternatives)))
         print(f"{line}\\t{alternatives}", file=outfile)

   """
}


/*
 * Extract cell barcodes and UMIs
 */
process UMITOOLS_EXTRACT {

   tag "${sample_id}"

   errorStrategy 'retry'
   maxRetries 2

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.extract.log" )

   input:
   tuple val( sample_id ), path( reads ), val( umis ), path( whitelist )

   output:
   tuple val( sample_id ), path( "*.gz" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   umi_tools extract \
		--bc-pattern "${umis[0]}" \
		--bc-pattern2 "${umis[1]}" \
      --extract-method regex \
      --error-correct-cell \
      --quality-encoding phred33 \
      --whitelist ${whitelist} \
      --log ${sample_id}.extract.log \
		--stdin ${reads[0]} \
		--read2-in ${reads[1]} \
		--stdout ${sample_id}_R1.extracted.fastq.gz \
		--read2-out ${sample_id}_R2.extracted.fastq.gz
   """
}

/*
 * Index the reference genome for use by Bowtie2.
 */
process BOWTIE2_INDEX {

   tag "${genome_id}"
   
   input:
   tuple val( genome_id ), path( fasta )

   output:
   tuple val( genome_id ), path( "*.bt2" )

   script:
   """
   bowtie2-build ${fasta} ${genome_id}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BOWTIE2_ALIGN {

   tag "${sample_id}-${genome_acc}" 

   errorStrategy 'retry'
   maxRetries 2

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads ), path( idx ), val( genome_acc )

   output:
   tuple val( sample_id ), path( "*.bam" ), emit: main
   path "*.log", emit: logs

   script:
   """
   bowtie2 \
      -x ${genome_acc} \
      --rdg 10,10 \
      --very-sensitive-local \
      -1 ${reads[0]} -2 ${reads[1]} \
      -S ${sample_id}.mapped.sam \
      2> ${sample_id}.bowtie2.log
   samtools sort ${sample_id}.mapped.sam \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   rm ${sample_id}.mapped.sam
   """
}

/*
 * Dedupicate reads based on mapping coordinates and UMIs
 */
process UMITOOLS_DEDUP {

   tag "${sample_id}" 

   label 'big_mem'
   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( processed_o, 
               mode: 'copy' )
   input:
   tuple val( sample_id ), path( bamfile )

   output:
   tuple val( sample_id ), path( "*.bam" ), emit: main
   path( "*.log" ), emit: logs
   path( "*.tsv"), emit: stats

   script:
   """
   samtools index ${bamfile}
   umi_tools dedup \
		--per-cell \
      --paired \
      --method unique \
      --chimeric-pairs discard \
      --unmapped-reads discard \
		--stdin ${bamfile} \
		--log ${sample_id}.dedup.log \
      --stdout ${sample_id}.dedup.bam
   #umicollapse bam -i ${bamfile} -o ${sample_id}.dedup.bam --paired --two-pass
   """
}


process UMICOLLAPSE {

   tag "${sample_id}" 

   label 'big_mem'
   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( processed_o, 
               mode: 'copy' )
   input:
   tuple val( sample_id ), path( bamfile ), path( umicollapse_repo )

   output:
   tuple val( sample_id ), path( "*.bam" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   samtools index ${bamfile}
   java -jar ${umicollapse_repo}/umicollapse.jar bam -i ${bamfile} -o ${sample_id}.dedup.bam --paired --two-pass 2>&1 > ${sample_id}.umicollapse.log
   """
}


/*
 * Use featureCounts (from the SubReads package) to count transcripts mapping to each gene
 * ## $(ANN_TYPE): which of "CDS", "gene", "mRNA", etc
 * ## $(LABEL): the tag from column 9 to use to label transcript counts, e.g. "Locus", "Name"
 */
process FEATURECOUNTS {

   tag "${sample_id}"

   errorStrategy 'retry'
   maxRetries 2

   publishDir( counts_o, 
               mode: 'copy', 
               pattern: "*.tsv" )
   publishDir( counts_o, 
               mode: 'copy', 
               pattern: "*.summary" )

   input:
   tuple val( sample_id ), path( bamfile ), path( gff )

   output:
   tuple val( sample_id ), path( "*.bam" ), emit: main
   tuple path( "*.summary" ), path( "*.tsv" ), emit: logs

   script:
   """
   samtools index ${bamfile}
   featureCounts \
      -p -C -s ${params.strand} \
      -d ${params.min_length} -D 10000  \
      -t ${params.ann_type} -g ${params.label} \
      -a ${gff} \
      -R BAM \
      -T ${task.cpus} \
      --countReadPairs \
      --verbose \
      --extraAttributes gene_biotype,locus_tag \
      -o ${sample_id}.featureCounts.tsv ${bamfile}
   """
}

/*
 * Count unique UMIs per cell per gene
 */
process UMITOOLS_COUNT {

   tag "${sample_id}"

   errorStrategy 'retry'
   maxRetries 2

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile )

   output:
   tuple val( sample_id ), path( "*.tsv" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   samtools sort ${bamfile} -o ${bamfile.getBaseName()}.sorted.bam
   samtools index ${bamfile.getBaseName()}.sorted.bam
   umi_tools count \
		--per-gene \
      --per-cell \
		--gene-tag XT \
      --log ${sample_id}.count.log \
		--stdin ${bamfile.getBaseName()}.sorted.bam \
      --stdout ${sample_id}.umitools_count.tsv
   rm ${bamfile.getBaseName()}.sorted.bam
   """
}

/*
 * Make log report
 */
process MULTIQC {

   errorStrategy 'retry'
   maxRetries 2

   publishDir( multiqc_o, 
               mode: 'copy' )

   input:
   path( '*' )

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/