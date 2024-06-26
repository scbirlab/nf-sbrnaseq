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

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println """\
         S C B I R   s b R N A - S E Q   P I P E L I N E
         ===============================================
         Nextflow pipeline to process demultiplexed Illumina paired-end 
         FASTQ files from multiple bacterial samples into a gene x
         cell count table.

         Usage:
            nextflow run sbcirlab/nf-sbrnaseq --sample_sheet <csv> --fastq_dir <dir> --genome_fasta_dir <dir> --genome_gff_dir <dir>
            nextflow run sbcirlab/nf-sbrnaseq -c <config-file>

         Required parameters:
            sample_sheet         Path to a CSV containing sample IDs matched with FASTQ filenames, genome information, and adapter sequences.
            fastq_dir            Path to directory containing the FASTQ file.
            genome_fasta_dir     Path to directory containing genome FASTA files (for mapping)
            genome_gff_dir       Path to directory containing genome GFF files (for feature counting)

         Optional parameters (with defaults):   
            trim_qual = 10       For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = 11      For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
            umitools_error = 6   For `umitools`, the number of errors allowed to correct cell barcodes
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
if ( !params.genome_fasta_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to genome_fasta_dir")
}
if ( !params.genome_gff_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to genome_gff_dir")
}

wd = file(params.sample_sheet)
working_dir = wd.getParent()

processed_o = "${working_dir}/outputs/processed"
counts_o = "${working_dir}/outputs/counts"
multiqc_o = "${working_dir}/outputs/multi_qc"

log.info """\
         S C B I R   s b R N A - S E Q   P I P E L I N E
         ===============================================
         inputs
            sample sheet   : ${params.sample_sheet}
            fastq directory: ${params.fastq_dir}
         genome locations
            FASTA          : ${params.genome_fasta_dir}
            GFF            : ${params.genome_gff_dir}
         trimming 
            quality        : ${params.trim_qual}
            minimum length : ${params.min_length}
         UMI detection
            Error number   : ${params.umitools_error}
         output            
            Processed      : ${processed_o}
            Counts         : ${counts_o}
            MultiQC        : ${multiqc_o}
         """
         .stripIndent()

dirs_to_make = [processed_o, counts_o, multiqc_o]

log.info  """
         Making directories: 
          """.stripIndent()

dirs_to_make.each { 
   log.info "${it}: " 
   log.info file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}

/*
========================================================================================
   Create Channels
========================================================================================
*/

csv_ch = Channel.fromPath( params.sample_sheet, 
                           checkIfExists: true )
               .splitCsv( header: true )

sample_ch = csv_ch.map { tuple( it.genome_id,
                                it.sample_id, 
                                tuple( it.adapter_read1_3prime,
                                       it.adapter_read2_3prime,
                                       it.adapter_read1_5prime,
                                       it.adapter_read2_5prime ),
                                tuple( it.umi_read1, 
                                       it.umi_read2 ),
                                file( "${params.fastq_dir}/*${it.fastq_pattern}*",
                                      checkIfExists: true ).sort() ) }

genome_ch = csv_ch.map { tuple( it.genome_id,
                                file( "${params.genome_fasta_dir}/${it.genome_id}.fasta",
                                      checkIfExists: true),
                                file( "${params.genome_gff_dir}/${it.genome_id}.gff3",
                                      checkIfExists: true) ) }
                  .unique()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   sample_ch | FASTQC 

   sample_ch | TRIM_CUTADAPT 
   TRIM_CUTADAPT.out.main | UMITOOLS_WHITELIST
   // TRIM_CUTADAPT.out.main.join( UMITOOLS_WHITELIST.out.main, 
   //                              by: 1 ) | UMITOOLS_EXTRACT
   TRIM_CUTADAPT.out.main | UMITOOLS_EXTRACT

   genome_ch | BOWTIE2_INDEX 

   UMITOOLS_EXTRACT.out.main.combine( BOWTIE2_INDEX.out, 
                                      by: 0 ) | BOWTIE2_ALIGN

   BOWTIE2_ALIGN.out.main \
      | UMITOOLS_DEDUP

   UMITOOLS_DEDUP.out.main \
      | FEATURECOUNTS

   FEATURECOUNTS.out.main \
      | UMITOOLS_COUNT 

   TRIM_CUTADAPT.out.logs.concat(
         FASTQC.out.logs,
         BOWTIE2_ALIGN.out.logs,
         FEATURECOUNTS.out.logs )
      .flatten()
      .unique()
      .collect() \
      | MULTIQC

}

/* 
 * Do quality control checks
 */
process FASTQC {

   memory '32GB'

   tag "${sample_id}" 

   input:
   tuple val( genome_id ), val( sample_id ), val( adapters ), val( umis ), path( 'reads_*' )

   output:
   path( "*.zip" ), emit: logs

   script:
   """
   zcat reads_* > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads 4 ${sample_id}.fastq
   rm ${sample_id}.fastq
   """

}

/*
 * Trim adapters from reads
 */
process TRIM_CUTADAPT {

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( genome_id ), val( sample_id ), val( adapters ), val( umis ), path( reads, stageAs: "?/*" )

   output:
   tuple val( genome_id ), val( sample_id ), val( adapters ), val( umis ), path( "*.with-adapters.fastq.gz" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   for i in \$(seq 1 2)
   do
      zcat */*_R"\$i"*.fastq.gz | gzip --best > ${sample_id}_polyA_R"\$i".fastq.gz
   done

   cutadapt \
		-a "A{10}" \
      -A "T{10}" \
		-q ${params.trim_qual} \
		--nextseq-trim ${params.trim_qual} \
		--minimum-length ${params.min_length} \
		--report full \
      --action trim \
		-o ${sample_id}_3prime_R1.fastq.gz \
      -p ${sample_id}_3prime_R2.fastq.gz \
		${sample_id}_polyA_R?.fastq.gz > ${sample_id}.polyA.cutadapt.log

   cutadapt \
		-a "${adapters[0]}" \
      -A "${adapters[1]}" \
      --minimum-length ${params.min_length} \
		--report full \
      --action trim \
		-o ${sample_id}_5prime_R1.fastq.gz \
      -p ${sample_id}_5prime_R2.fastq.gz \
		${sample_id}_3prime_R?.fastq.gz > ${sample_id}.3prime.cutadapt.log

   cutadapt \
		-g "${adapters[2]}" \
      -G "${adapters[3]}" \
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

/*
 * Extract cell barcodes and UMIs
 */
process UMITOOLS_EXTRACT {

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.extract.log" )

   input:
   tuple val( genome_id ), val( sample_id ),  val( adapters ), val( umis ), path( reads )

   output:
   tuple val( genome_id ), val( sample_id ), path( "*.gz" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   umi_tools extract \
		--bc-pattern "${umis[0]}" \
		--bc-pattern2 "${umis[1]}" \
      --extract-method regex \
      --error-correct-cell \
      --quality-filter-mask ${params.trim_qual} \
      --quality-encoding phred33 \
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
   tuple val( genome_id ), path( fasta ), path( gff )

   output:
   tuple val( genome_id ), path( "*.bt2" ), path( gff )

   script:
   """
   bowtie2-build ${fasta} ${genome_id}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BOWTIE2_ALIGN {

   tag "${sample_id}-${genome_id}" 

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( genome_id ), val( sample_id ), path( reads ), path( idx ), path( gff )

   output:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   bowtie2 \
      -x ${genome_id} \
      --rdg 10,10 \
      --very-sensitive-local \
      -1 ${reads[0]} -2 ${reads[1]} \
      -S ${sample_id}.mapped.sam \
      2> ${sample_id}.bowtie2.log
   samtools sort ${sample_id}.mapped.sam \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   samtools index ${sample_id}.mapped.bam
   rm ${sample_id}.mapped.sam
   """
}

/*
 * Dedupicate reads based on mapping coordinates and UMIs
 */
process UMITOOLS_DEDUP {
   errorStrategy 'ignore'

   tag "${sample_id}" 

   publishDir( processed_o, 
               mode: 'copy')
   input:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( bamfile ), path( bam_idx )

   output:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( "*.bam" ), emit: main
   path( "*.log" ), emit: logs
   path( "*.tsv"), emit: stats

   script:
   """
   umi_tools dedup \
		--per-cell \
      --output-stats ${sample_id}.dedup \
		--stdin ${bamfile} \
		--log ${sample_id}.dedup.log \
      --stdout ${sample_id}.dedup.bam
   """
}


/*
 * Use featureCounts (from the SubReads package) to count transcripts mapping to each gene
 * ## $(ANN_TYPE): which of "CDS", "gene", "mRNA", etc
 * ## $(LABEL): the tag from column 9 to use to label transcript counts, e.g. "Locus", "Name"
 */
process FEATURECOUNTS {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy', 
               pattern: "*.tsv" )
   publishDir( counts_o, 
               mode: 'copy', 
               pattern: "*.summary" )

   input:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( bamfile )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   tuple path( "*.summary" ), path( "*.tsv" ), emit: logs

   script:
   """
   featureCounts \
      -p -B -C -s ${params.strand} \
      -d ${params.min_length} -D 10000  \
      -t ${params.ann_type} -g ${params.label} \
      -a ${gff} \
      -R BAM \
      -o ${sample_id}.featureCounts.tsv ${bamfile}
   samtools index *featureCounts.bam
   """
}

/*
 * Count unique UMIs per cell per gene
 */
process UMITOOLS_COUNT {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( bam_idx )

   output:
   tuple val( sample_id ), path( "*.tsv" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   umi_tools count \
		--per-gene --per-cell \
		--gene-tag XT \
      --log ${sample_id}.count.log \
		--stdin  ${bamfile} \
      --stdout ${sample_id}.umitools_count.tsv
   """
}

/*
 * Make log report
 */
process MULTIQC {

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