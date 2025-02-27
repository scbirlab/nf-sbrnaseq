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
            nextflow run sbcirlab/nf-sbrnaseq --sample_sheet <csv> [--fastq-dir <dir>|--from-sra]
            nextflow run sbcirlab/nf-sbrnaseq -c <config-file>

         Required parameters:
            sample_sheet               Path to a CSV containing sample IDs matched with FASTQ filenames, genome information, and adapter sequences.

         If using local FASTQ data (the default behavior):
            fastq_dir                  Path to directory containing the FASTQ file.
          
         Optional parameters (with defaults):  
            from_sra = false           Whether to fetch FASTQ data from the SRA.
            allow_cell_errors = true   Whether to allow 1 error when matching cell barcodes in the whitelist.
            trim_qual = 5              For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = "9:38"        For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
            strand = 1                 For `featureCounts`, the strandedness of RNA-seq. `1` for forward, `2` for reverse.
            ann_type = 'gene'          For `featureCounts`, features from GFF column 3 to use for counting
            label = 'Name'             For `featureCounts`, one or more (comma-separated) fields from column 9 of GFF for labeling counts

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
if ( !params.from_sra ) {
   if ( !params.fastq_dir ) {
      throw new Exception("!!! PARAMETER MISSING: Please provide a path to fastq_dir")
   }
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
   if ( params.from_sra ) {
      Channel.of( params.ncbi_api_key ).set { ncbi_api_key }
      csv_ch
         .map { tuple( 
               it.sample_id,
               it.Run
            ) }
         .combine( ncbi_api_key ) 
         | PULL_FASTQ_FROM_SRA
      PULL_FASTQ_FROM_SRA.out
         .set { reads_ch }  // sample_id, reads
   } else {
      csv_ch
         .map { tuple( 
            it.sample_id,
            file( 
               "${params.fastq_dir}/*${it.fastq_pattern}*",
               checkIfExists: true 
            ).sort()
         ) }
         .set { reads_ch }  // sample_id, [reads]
   }
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
   if ( ! params.use_star ) {
      DOWNLOAD_GENOME.out
         .map { it[0..1] }  // genome_acc, genome
         | BOWTIE2_INDEX   // genome_acc, [genome_idx]
      BOWTIE2_INDEX.out.set { genome_idx0 }
   } else {
      DOWNLOAD_GENOME.out
         | STAR_INDEX   // genome_acc, [genome_idx]
      STAR_INDEX.out.set { genome_idx0 }
   }
   genome_ch
      .map { it[1..0] }  // genome_acc, sample_id
      .combine( genome_idx0, by: 0 )  // genome_acc, sample_id, [genome_idx]
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
   
   if ( params.allow_cell_errors ) {
      BUILD_WHITELIST.out |
       ADD_WHITELIST_ERRORS  // wl_id, whitelist_err
      ADD_WHITELIST_ERRORS.out
       .set { whitelist0 }
   } else {
      BUILD_WHITELIST.out.set { whitelist0 }
   }
   barcode_ch
      .map { it[1..0] }  // wl_id, sample_id
      .combine( whitelist0, by: 0 )  // wl_id, sample_id, whitelist_err
      .map { it[1..-1] }  // sample_id, whitelist_err
      .set { whitelists }

   TRIM_CUTADAPT.out.main  // sample_id, [reads]
      .combine( umi_ch, by: 0 )  // sample_id, [reads], umis
      .combine( whitelists, by: 0 )  // sample_id, [reads], umis, whitelist
      | UMITOOLS_EXTRACT  // sample_id, [reads]

   UMITOOLS_EXTRACT.out.main
      .combine( genome_idx, by: 0 )  // sample_id, [reads], [genome_idx], genome_acc
      | ( params.use_star ? STAR_ALIGN : BOWTIE2_ALIGN ) // sample_id, [bam_bai]
      | set { mapped_reads }

   mapped_reads.main
      .combine( DOWNLOAD_UMICOLLAPSE.out )  // sample_id, [bam_bai], umicollapse_repo
      | UMICOLLAPSE  // sample_id, dedup_bam

   UMICOLLAPSE.out.main
      .combine( genome_gff, by: 0 )  // sample_id, dedup_bam, gff
      | FEATURECOUNTS  // sample_id, [count_bam_bai]

   FEATURECOUNTS.out.main
      | UMITOOLS_COUNT   // sample_id, counts
   // UMITOOLS_COUNT.out.main
   //    | MAKE_COUNT_MATRIX
   UMITOOLS_COUNT.out.main
      .combine( FEATURECOUNTS.out.table, by: 0 )  // sample_id, umitools_counts, featurecounts_counts
      | JOIN_FEATURECOUNTS_UMITOOLS
   JOIN_FEATURECOUNTS_UMITOOLS.out.all_counts
      | PLOT_UMI_DISTRIBUTIONS   
   JOIN_FEATURECOUNTS_UMITOOLS.out.cells_per_gene
      | PLOT_CELLS_PER_GENE_DISTRIBUTION
   JOIN_FEATURECOUNTS_UMITOOLS.out.genes_per_cell
      | PLOT_GENES_PER_CELL_DISTRIBUTION

   JOIN_FEATURECOUNTS_UMITOOLS.out.all_counts
      | MAKE_ANNDATA

   TRIM_CUTADAPT.out.logs
      .concat(
         FASTQC.out.multiqc_logs,
         mapped_reads.logs,
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
   tuple val( accession ), path( "all-nucleotides.fna" ), path( "all-annotations.gff" )

   script:
   """
   ACCESSIONS=\$(echo "${accession}" | tr '+' ' ')
   echo "\$ACCESSIONS"
   WEB_ROOT="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
   WEB_TAIL="download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
   for acc in \$ACCESSIONS
   do
      wget "\$WEB_ROOT/\$acc/\$WEB_TAIL" \
         -O \$acc"_genome-out"
      unzip -o \$acc"_genome-out" "ncbi_dataset/data/\$acc/"{"\$acc"_*_genomic.fna,*.gff}
      mv ncbi_dataset/data/*/"\$acc"_*_genomic.fna \$acc.fna
      mv ncbi_dataset/data/*/*.gff \$acc.gff
   done

   cat *.fna > all-nucleotides.fna
   cat *.gff > all-annotations.gff
   """
}


// Get FASTQ
process PULL_FASTQ_FROM_SRA {

   tag "${sample_id}-${sra_run_id}" 

   label 'big_mem'
   time '24 h'

   input:
   tuple val( sample_id ), val( sra_run_id ), val( ncbi_api_key )

   output:
   tuple val( sample_id ), path( "*.with-idx_R?.fastq.gz" )

   script:
   """
   NCBI_API_KEY=${ncbi_api_key} \
   fastq-dump \
      --read-filter pass \
      --origfmt --defline-seq '@rd.\$si:\$sg:\$sn' \
      --split-3 ${sra_run_id}

   for i in \$(seq 1 2)
   do
      if [ \$i -eq 1 ]
      then
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4==1 { a=\$2; alen=length(a); print \$0 } NR%4==2 { print a \$0 } NR%4==3 { print "+" } NR%4==0 { s = sprintf("%*s", alen, ""); print gensub(".", "F", "g", s) \$0 }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      else
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4<3 { print \$0 } NR%4==3 { print "+" }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      fi
   done
   rm *.fastq
   """

   stub:
   """
   NCBI_API_KEY=${ncbi_api_key} \
   fastq-dump \
      -X 1000000 \
      --read-filter pass \
      --origfmt --defline-seq '@rd.\$si:\$sg:\$sn' \
      --split-3 ${sra_run_id}

   for i in \$(seq 1 2)
   do
      if [ \$i -eq 1 ]
      then
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4==1 { a=\$2; alen=length(a); print \$0 } NR%4==2 { print a \$0 } NR%4==3 { print "+" } NR%4==0 { s = sprintf("%*s", alen, ""); print gensub(".", "F", "g", s) \$0 }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      else
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4<3 { print \$0 } NR%4==3 { print "+" }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      fi
   done
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
   label "big_mem" 
   
   // errorStrategy 'retry'
   // maxRetries 1

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads, stageAs: "?/*" ), val( adapters5 ), val( adapters3 )
   tuple val( trim_qual ), val( min_length )

   output:
   tuple val( sample_id ), path( "*.with-adapters.fastq.gz" ), emit: main
   path( "*.log" ), emit: logs

   script:
   def adapter_3prime_R = (adapters3[1].length() > 0 ? "-A '${adapters3[1]}'" : "")
   def adapter_5prime_R = (adapters5[1].length() > 0 ? "-G '${adapters5[1]}'" : "")
   def adapter_5prime_F = (adapters5[0].length() > 0 ? "-g '${adapters5[0]}'" : "")
   """
   for i in \$(seq 1 2)
   do
      cat */*_R"\$i"*.fastq.gz > ${sample_id}_polyA_R"\$i".fastq.gz
   done

   cutadapt \
		-a "A{10}" \
		-q ${trim_qual} \
		--nextseq-trim ${trim_qual} \
		--minimum-length ${min_length} \
		--report full \
      --action trim \
      -j ${task.cpus} \
		-o ${sample_id}_3prime_R1.fastq.gz \
      -p ${sample_id}_3prime_R2.fastq.gz \
		${sample_id}_polyA_R?.fastq.gz > ${sample_id}.polyA.cutadapt.log

   cutadapt \
		-a "${adapters3[0]}" \
      ${adapter_3prime_R} \
      --minimum-length ${min_length} \
		--report full \
      --action trim \
      -j ${task.cpus} \
		-o ${sample_id}_5prime_R1.fastq.gz \
      -p ${sample_id}_5prime_R2.fastq.gz \
		${sample_id}_3prime_R?.fastq.gz > ${sample_id}.3prime.cutadapt.log

   ADAPT5_ALL="${adapter_5prime_F}${adapter_5prime_R}"
   ADAPT5_LEN=\${#ADAPT5_ALL}
   if [ \$ADAPT5_LEN -gt 0 ]
   then
      cutadapt \
         ${adapter_5prime_F} ${adapter_5prime_R} \
         --report full \
         --action retain \
         --discard-untrimmed \
         -j ${task.cpus} \
         -o ${sample_id}_R1.with-adapters.fastq.gz \
         -p ${sample_id}_R2.with-adapters.fastq.gz \
         ${sample_id}_5prime_R?.fastq.gz > ${sample_id}.5prime.cutadapt.log
   else
      for i in \$(seq 1 2)
      do
         mv ${sample_id}_5prime_R\$i.fastq.gz ${sample_id}_R\$i.with-adapters.fastq.gz
      done
   fi

   rm ${sample_id}_{3prime,polyA}_R?.fastq.gz
   """
}

// Build whitelist from provided barcode files
process BUILD_WHITELIST {
   tag "${whitelist_id}"

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
   label "big_time"

   input:
   tuple val( whitelist_id ), path( whitelist )

   output:
   tuple val( whitelist_id ), path( "*.whitelist-err.txt" )

   script:
   """
   #!/usr/bin/env python

   from itertools import product

   ALPHABET = "ATCGN"
   INPUT_FILE = "${whitelist}"
   OUTPUT_FILE = "${whitelist_id}.whitelist-err.txt"

   with open(OUTPUT_FILE, 'w') as outfile, open(INPUT_FILE, 'r') as infile:
      for line in infile:
         bc_length = len(line)
         break
      infile.seek(0)
      bc_too_long = (bc_length * (len(ALPHABET) - 1)) > 120
      if not bc_too_long:
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
      else:
         for line in infile:
            print(line.strip(), file=outfile)

   """
}


/*
 * Extract cell barcodes and UMIs
 */
process UMITOOLS_EXTRACT {

   tag "${sample_id}"
   label "big_mem"
   time "2d"  // in case of very large whitelists

   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.extract.log" )

   input:
   tuple val( sample_id ), path( reads ), val( umis ), path( whitelist )

   output:
   tuple val( sample_id ), path( "*.extracted.fastq.gz" ), emit: main
   path "*.log", emit: logs

   script:
   def error_correct = ( params.allow_cell_errors ? "--error-correct-cell" : "" )
   """
   MAX_SIZE=500000000  # 500 million
   split -l \$MAX_SIZE "${whitelist}" chunk_

   for f in chunk_*
   do
      umi_tools extract \
         --bc-pattern "${umis[0]}" \
         --bc-pattern2 "${umis[1]}" \
         --extract-method regex ${error_correct} \
         --quality-encoding phred33 \
         --whitelist "\$f" \
         --log ${sample_id}.extract.log \
         --stdin ${reads[0]} \
         --read2-in ${reads[1]} \
         --stdout "\$f"_R1.fastq.gz \
         --read2-out "\$f"_R2.fastq.gz
   done

   for i in \$(seq 1 2)
   do
      cat chunk_*_R\$i.fastq.gz > ${sample_id}_R\$i.extracted.fastq.gz
   done
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
   label "big_mem"

   errorStrategy 'retry'
   maxRetries 1

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
      -L 13 -N 1 \
      --very-sensitive-local \
      -p ${task.cpus} \
      -1 ${reads[0]} -2 ${reads[1]} \
      -S ${sample_id}.mapped.sam \
      2> ${sample_id}.bowtie2.log
   samtools sort -@ ${task.cpus} ${sample_id}.mapped.sam \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   rm ${sample_id}.mapped.sam
   """
}


/*
 * Index the reference genome for use by STAR.
 */
process STAR_INDEX {

   tag "${genome_id}"
   
   input:
   tuple val( genome_id ), path( fasta ), path( gff )

   output:
   tuple val( genome_id ), path( "star-index" )

   script:
   """
   mkdir star-index
   STAR \
      --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir star-index \
      --genomeFastaFiles ${fasta} \
      --sjdbOverhang 100 \
      --sjdbGTFfile ${gff} \
      --sjdbGTFfeatureExon CDS \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbGTFtagExonParentGene Parent \
      --sjdbGTFtagExonParentGeneName gene \
      --sjdbGTFtagExonParentGeneType gene_biotype
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process STAR_ALIGN {

   tag "${sample_id}-${genome_acc}" 
   label "big_mem"
   time "2d"

   // errorStrategy 'retry'
   // maxRetries 1

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads ), path( idx ), val( genome_acc )

   output:
   tuple val( sample_id ), path( "*.bam" ), emit: main
   path "*.log", emit: logs

   script:
   """
   # make sure no 0-length reads
   for f in ${reads}
   do
      zcat \$f \
      | awk 'NR%4 == 1 { name = \$0 } NR%4 == 2 { if (length(\$0) > 0) { seq = \$0 } else { seq = "N" } } NR%4 == 0 { if (length(\$0) > 0) { qual = \$0 } else { qual = "?" } } ( NR > 1 && NR%4 == 1 ) { print name; print seq; print "+"; print qual } END { print name; print seq; print "+"; print qual }' \
      | gzip --best \
      > \$(basename \$f .fastq.gz).padded.fastq.gz
   done

   STAR \
      --runThreadN ${task.cpus} \
      --genomeDir ${idx} \
      --readFilesIn *_R1*.padded.fastq.gz *_R2*.padded.fastq.gz \
      --readFilesCommand zcat \
      --alignEndsType Extend5pOfRead1 \
      --alignIntronMax 1 \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNmax 50 \
      --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
      --twopassMode Basic \
      --quantMode GeneCounts \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM ${task.memory.getBytes()} \
      --outBAMsortingThreadN ${task.cpus}
   mv Log.final.out ${sample_id}.star.log
   mv Aligned.sortedByCoord.out.bam ${sample_id}.mapped.bam
   rm *.padded.fastq.gz
   """
}


process UMICOLLAPSE {

   tag "${sample_id}" 

   label 'med_mem'
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
   java -jar -Xmx${task.memory.getMega()}m -Xss1024m ${umicollapse_repo}/umicollapse.jar bam \
      -i ${bamfile} -o ${sample_id}.dedup.bam --paired --two-pass 2>&1 > ${sample_id}.umicollapse.log
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
               pattern: "*.log" )

   input:
   tuple val( sample_id ), path( bamfile ), path( gff )

   output:
   tuple val( sample_id ), path( "*.bam" ), emit: main
   tuple val( sample_id ), path( "*.tsv" ), emit: table
   path "*.summary", emit: logs

   script:
   """
   samtools index ${bamfile}
   featureCounts \
      -p -C -s ${params.strand} \
      -t ${params.ann_type} -g ${params.label} \
      -a ${gff} \
      -R BAM \
      -T ${task.cpus} \
      --countReadPairs \
      --verbose \
      --extraAttributes Name,gene_biotype,locus_tag \
      -o ${sample_id}.featureCounts.tsv ${bamfile}
   cp ${sample_id}.featureCounts.tsv.summary ${sample_id}.featureCounts.log
   """
}

/*
 * Count unique UMIs per cell per gene
 */
process UMITOOLS_COUNT {

   tag "${sample_id}"

   errorStrategy 'retry'
   maxRetries 2

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
      --stdout ${sample_id}.umitools_count.tsv \
   && rm ${bamfile.getBaseName()}.sorted.bam
   """
}


process MAKE_COUNT_MATRIX {
   tag "${sample_id}"
   label "med_mem"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( umitools_table )

   output:
   tuple val( sample_id ), path( "*.tsv.gz" )

   script:
   """
   #!/usr/bin/env python
   from scipy.sparse import csr_matrix
   from pandas.api.types import CategoricalDtype
   import pandas as pd

   df = pd.read_csv("${umitools_table}", sep='\\t')
   cell_cat = CategoricalDtype(sorted(df["cell"].unique()), ordered=True)
   gene_cat = CategoricalDtype(sorted(df["gene"].unique()), ordered=True)

   cell_idx = df["cell"].astype(cell_cat).cat.codes
   gene_idx = df["gene"].astype(gene_cat).cat.codes
   matrix = csr_matrix(
      (df["count"], (cell_idx, gene_idx)),
      shape=(cell_cat.categories.size, gene_cat.categories.size)
   )
   df = pd.DataFrame.sparse.from_spmatrix(
      matrix,
      index=cell_cat.categories,
      columns=gene_cat.categories,
   )
   df.to_csv("${sample_id}.count-matrix.tsv.gz", sep='\\t')
   """
}


process MAKE_ANNDATA {
   tag "${sample_id}"
   label "med_mem"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( counts_table )

   output:
   tuple val( sample_id ), path( "*.h5ad" ), path( "*.png" ), path( "figures/*.png" ), emit: main
   path "*.log", emit: logs

   script:
   """
   ${projectDir}/bin/make_anndata.py ${sample_id} ${counts_table}
   """
}



process JOIN_FEATURECOUNTS_UMITOOLS {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( umitools_table ), path( featurecounts_table )

   output:
   tuple val( sample_id ), path( "*.all-counts.tsv" ), emit: all_counts
   tuple val( sample_id ), path( "*.chr-per-cell*.tsv" ), emit: chr_per_cell
   tuple val( sample_id ), path( "*.umis-and-cells-per-gene.tsv" ), emit: cells_per_gene
   tuple val( sample_id ), path( "*.genes-per-cell.tsv" ), path( "*.genes-per-biotype-per-cell.tsv" ), emit: genes_per_cell

   script:
   """
   set -x
   grep -v '^#' "${featurecounts_table}" | tr \$'\\t' , > fc0.csv
   LOCUS_TAG_COL=\$(head -n1 fc0.csv | tr , \$'\\n' | grep -n "Geneid" | cut -d: -f1)
   cat \
   <(head -n1 fc0.csv | awk -F, -v OFS=, '{ print "gene_id",\$0 }' | sed 's/,${sample_id}.*\\.dedup\\.bam\$/,featurecounts_count/' )\
   <(tail -n+2 fc0.csv \
      | awk -F, -v OFS=, '{ print \$'\$LOCUS_TAG_COL',\$0}' \
      | sort -k1 -t, ) \
   > fc.csv

   cat ${umitools_table} | tr \$'\\t' , > ut0.csv
   cat \
   <(head -n1 ut0.csv | sed 's/,count\$/,umi_count/;s/^gene,/gene_id,/;s/,cell,/,cell_barcode,/' ) \
   <(tail -n+2 ut0.csv | sort -k1 -t, ) \
   > ut.csv

   python -c \
      'import pandas as pd; pd.read_csv("fc.csv").merge(pd.read_csv("ut.csv")).assign(sample_id="${sample_id}").set_index(["sample_id", "gene_id", "locus_tag", "cell_barcode"]).reset_index().to_csv("${sample_id}.all-counts.tsv", sep="\\t", index=False)'
   
   # count the number of chromosomes per cell
   printf 'n_chromosomes\\tcell_barcode\\n' >  "${sample_id}.chr-per-cell.tsv"
   cat "${sample_id}.all-counts.tsv" \
   | tail -n+2 \
   | cut -f4,6 | sort -k1 \
   | uniq | cut -f1 \
   | uniq -c | sort -k1 -n \
   >> "${sample_id}.chr-per-cell.tsv"
   printf 'n_chromosomes\\tn_cell_barcodes\\n' >  "${sample_id}.chr-per-cell-summary.tsv"
   awk 'NR > 1 { print \$1 }' "${sample_id}.chr-per-cell.tsv" \
   | uniq -c \
   >> "${sample_id}.chr-per-cell-summary.tsv"

   python -c \
      'import pandas as pd; df = pd.read_csv("ut.csv").groupby("gene_id"); df[["cell_barcode"]].count().reset_index().merge(df[["umi_count"]].sum().reset_index()).merge(pd.read_csv("fc.csv")).rename(columns=dict(cell_barcode="cell_count")).to_csv("${sample_id}.umis-and-cells-per-gene.tsv", sep="\\t", index=False)'
   python -c \
      'import pandas as pd; df = pd.read_csv("${sample_id}.all-counts.tsv", sep="\\t", low_memory=False).groupby("cell_barcode"); df[["gene_id"]].count().reset_index().merge(df[["umi_count"]].sum().reset_index()).rename(columns=dict(gene_id="gene_count")).sort_values("umi_count").to_csv("${sample_id}.genes-per-cell.tsv", sep="\\t", index=False)' 
   python -c \
      'import pandas as pd; df = pd.read_csv("${sample_id}.all-counts.tsv", sep="\\t", low_memory=False).groupby(["cell_barcode", "gene_biotype"]); df[["gene_id"]].count().reset_index().merge(df[["umi_count"]].sum().reset_index()).rename(columns=dict(gene_id="gene_count")).sort_values("umi_count").to_csv("${sample_id}.genes-per-biotype-per-cell.tsv", sep="\\t", index=False)'
   """
}


process PLOT_UMI_DISTRIBUTIONS {
   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( umi_table )

   output:
   tuple val( sample_id ), path( "*.png" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   biotype_blocklist = [
      #"SRP_RNA",
      #"antisense_RNA",
   ]
   df = pd.read_csv(
      "${umi_table}", 
      sep='\\t',
      low_memory=False,
   ).query("~gene_biotype.isin(@biotype_blocklist)")  # very low representation, breaks histogram bins

   fig, axes = scattergrid(
      df,
      grid_columns=["umi_count", "umi_count"],
      log=["umi_count", "umi_count"],
      grouping=["gene_biotype"],
      aspect_ratio=1.25,
   )
   for ax in fig.axes:
      ax.set(yscale="log")
   figsaver(format="png")(
      fig=fig,
      name='${sample_id}.umi-hist',
   )

   """
   
}


process PLOT_CELLS_PER_GENE_DISTRIBUTION {
   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( cells_per_gene_table )

   output:
   tuple val( sample_id ), path( "*.png" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   df = pd.read_csv(
      "${cells_per_gene_table}", 
      sep='\\t',
   ) 

   fig, axes = scattergrid(
      df,
      grid_columns=["umi_count", "cell_count"],
      log=["umi_count", "cell_count"],
   )
   figsaver(format="png")(
      fig=fig,
      name='${sample_id}.umis-and-cells-per-gene',
   )

   """
   
}


process PLOT_GENES_PER_CELL_DISTRIBUTION {
   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( gene_table ), path( gene_biotype_table )

   output:
   tuple val( sample_id ), path( "*.png" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   for f in ("${gene_table}", "${gene_biotype_table}"):
      df = pd.read_csv(
         f, 
         sep='\\t',
      ) 
      using_biotype = ("gene_biotype" in df)
      fig, axes = scattergrid(
         df,
         grid_columns=["umi_count", "gene_count"],
         log=["umi_count", "gene_count"],
         grouping=["gene_biotype"] if using_biotype else None,
         aspect_ratio=1.25 if using_biotype else 1.,
      )
      for ax in fig.axes:
         ax.set(yscale="log")
      figsaver(format="png")(
         fig=fig,
         name=f.split(".tsv")[0],
      )

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