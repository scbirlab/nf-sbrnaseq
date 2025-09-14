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

   Nextflow pipeline to process demultiplexed Illumina paired-end 
   FASTQ files from multiple bacterial samples into a gene x
   cell count table.
   """
   .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println pipeline_title + """\
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

log.info pipeline_title + """\
   Nanopore mode     : ${params.nanopore}
   inputs
      input dir.     : ${params.inputs}
      FASTQ dir.     : ${params.fastq_dir}
      sample sheet   : ${params.sample_sheet}
   trimming 
      quality        : ${params.trim_qual}
      minimum length : ${params.min_length}
   UMI-tools
      allow errors   : ${params.allow_cell_errors}
   Aligner           : ${params.mapper}
   FeatureCounts
      Strand         : ${params.strand}
      Annotation     : ${params.ann_type}
      Label          : ${params.label}
   UMIcollapse
      Source         : ${params.umicollapse_repo}
   output            : ${params.outputs}
   """
   .stripIndent()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

include { 
   bowtie2_index; 
   bowtie2_align;
} from './modules/bowtie.nf'
include { 
   join_featurecounts_UMItools;
   count_genomes_per_cell;
} from './modules/count_utils.nf'
include { featurecounts } from './modules/featurecounts.nf'
include { multiQC } from './modules/multiqc.nf'
include { 
   minimap_index; 
   minimap_align;
} from './modules/minimap.nf'
include { 
   fetch_genome_from_NCBI; 
   prefetch_from_SRA;
   download_FASTQ_from_SRA;
   prepend_reads_with_barcodes;
} from './modules/ncbi.nf'
include { 
   plot_UMI_distributions;
   plot_counts_per_gene; 
   plot_counts_per_cell;
} from './modules/plots.nf'
include { 
   gene_body_coverage;  
   gff2bed; 
} from './modules/rseqc.nf'
include { 
   STAR_index; 
   STAR_align;
} from './modules/star.nf'
include { 
   trim_using_cutadapt; 
   trim_nanopore_using_cutadapt; 
} from './modules/trimming.nf'
include { 
   fetch_UMIcollapse; 
   UMIcollapse;
} from './modules/umicollapse.nf'
include { 
   UMItools_count;
   UMItools_extract;
   concat_extractions;
   bam2table;
   counts_per_gene;
   counts_per_cell;
   counts_per_gene_per_cell;
} from './modules/umitools.nf'
include { 
   remove_multimappers;
   plot_bamstats;
   SAMtools_stats;
   SAMtools_coverage;
   SAMtools_flagstat;
 } from './modules/samtools.nf'
include { 
   build_AnnData;
   filter_AnnData;
   cluster_cells;
 } from './modules/scanpy.nf'
 include { 
   stringtie;
   merge_stringtie;
   stringtie_count;
 } from './modules/stringtie.nf'
include { fastQC } from './modules/qc.nf'
include { 
   add_errors_to_whitelist; 
   build_whitelist;
   chunk_whitelist;
} from './modules/whitelisting.nf'

workflow {

   /*
   ========================================================================================
      Read inputs
   ========================================================================================
   */

   Channel.fromPath( 
         params.sample_sheet, 
         checkIfExists: true 
      )
      .splitCsv( header: true )
      .set { csv_ch }

   if ( params.nanopore ) {

      csv_ch
         .map { tuple( 
            it.sample_id,
            it.adapter_5prime, 
            it.adapter_3prime,
         ) }
         .set { adapter_ch }  // sample_name, adapt5, adapt3

      csv_ch
         .map { tuple( 
            it.sample_id,
            tuple( it.umi )
         ) }
         .set { umi_ch }  // sample_id, umi

   } 
   
   else {

      csv_ch
         .map { tuple( 
            it.sample_id,
            tuple( it.adapter_read1_5prime, it.adapter_read2_5prime ),
            tuple( it.adapter_read1_3prime, it.adapter_read2_3prime ),
         ) }
         .unique()
         .set { adapter_ch }  // sample_name, [adapt5], [adapt3] 

      csv_ch
         .map { tuple( 
            it.sample_id,
            tuple( it.umi_read1, it.umi_read2 )
         ) }
         .unique()
         .set { umi_ch }  // sample_id, [umis]

   }

   csv_ch
      .map { tuple( 
         it.sample_id,
         tuple( it.bc1, it.bc2, it.bc3 )
            .collect {
               file( 
                  "${params.inputs}/${it}",
                  checkIfExists: true,
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
      csv_ch
         .map { tuple( it.sample_id, it.Run ) }
         | prefetch_from_SRA
         | download_FASTQ_from_SRA
         | prepend_reads_with_barcodes
      prepend_reads_with_barcodes.out
         .transpose()
         .groupTuple( by: 0 )
         .set { reads_ch }  // sample_id, [reads]
   } 
   
   else {
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
      .map { tuple( it.sample_id, it.genome_accession ) }
      .set { genome_ch }  // sample_id, genome_acc

   /*
   ========================================================================================
      Processing
   ========================================================================================
   */

   Channel.of( params.umicollapse_repo ) | fetch_UMIcollapse
   reads_ch | fastQC 
   genome_ch
      .map { it[1] }  // genome_acc
      .unique()
      | fetch_genome_from_NCBI   // genome_acc, genome, gff

   genome_ch
      .map { it[1..0] }  // genome_acc, sample_id
      .combine( fetch_genome_from_NCBI.out, by: 0 )  // genome_acc, sample_id, genome, gff
      .map { tuple( it[1], it[-1] ) }  // sample_id, gff
      .unique()
      .set { genome_gff }

   if ( !params.nanopore ) {

      trim_using_cutadapt(
         reads_ch.combine( adapter_ch, by: 0 ),  // sample_id, [reads], [adapt5], [adapt3]
         Channel.value( tuple( params.trim_qual, params.min_length ) )
      )  // sample_id, [reads]
      trim_using_cutadapt.out.main.set { trimmed }
      trim_using_cutadapt.out.logs.set { trim_logs }

   }

   else {

      trim_nanopore_using_cutadapt(
         reads_ch.combine( adapter_ch, by: 0 ),  // sample_id, [reads], adapt5, adapt3
         Channel.value( tuple( params.trim_qual, params.min_length ) )
      )  // sample_id, [reads]
      trim_nanopore_using_cutadapt.out.main.set { trimmed }
      trim_nanopore_using_cutadapt.out.logs.set { trim_logs }

   }
   
   build_whitelist(
      barcode_ch
         .map { it[1..2] }  // wl_id, [BCs]
         .unique(),
      Channel.value( params.reverse ),
   )
   
   if ( params.allow_cell_errors ) {
      
      build_whitelist.out 
         | add_errors_to_whitelist  // wl_id, whitelist_err
         | set { whitelist0 }

   } 
   
   else {
   
      build_whitelist.out
         | set { whitelist0 }
   
   }

   barcode_ch
      .map { it[1..0] }  // wl_id, sample_id
      .combine( 
         chunk_whitelist(
            whitelist0,
            Channel.value( params.whitelist_chunk_size ),
         ).transpose(), 
         by: 0,
       )  // wl_id, sample_id, whitelist_err
      .map { it[1..-1] }  // sample_id, whitelist_err
      .unique()
      .set { whitelists }

   trimmed  // sample_id, [reads]
      .combine( umi_ch, by: 0 )  // sample_id, [reads], umis
      .combine( whitelists, by: 0 )  // sample_id, [reads], umis, whitelist
      | UMItools_extract  // sample_id, [reads]
   UMItools_extract.out.main
      .transpose()
      .groupTuple( by: 0 )
      .map { tuple( it[0], it[1].sort() ) }
      .unique()
      | concat_extractions
      | set { extracted }

   if ( params.mapper == "star" ) {
      fetch_genome_from_NCBI.out
         | STAR_index 
         | set { genome_idx0 }
   }

   else if ( params.mapper == "bowtie2" ) {

      fetch_genome_from_NCBI.out
         .map { it[0..1] } 
         .unique()
         | bowtie2_index
         | set { genome_idx0 }

   } 

   else if ( params.mapper == "minimap2" ) {

      minimap_index(
         fetch_genome_from_NCBI.out
            .map { it[0..1] }
            .unique(),
         Channel.value( params.nanopore ),

      )
         | set { genome_idx0 }

   }

   else {
      error "Unsupported mapper: ${params.mapper}. Choose from: bowtie2, minimap2, star."
   }

   genome_idx0
      .combine( 
         genome_ch
            .map { it[1..0] }
            .unique(), 
         by: 0,
      )  // genome_acc, [genome_idx], sample_id
      .map { it[-1..0] }  // sample_id, [genome_idx], genome_acc
      .unique()
      .set { genome_idx }

   extracted
      .combine( genome_idx, by: 0 )  // sample_id, [reads], [genome_idx], genome_acc
      .set { pre_mapper }

   if ( params.mapper == "star" ) {

      pre_mapper
         | STAR_align
         | set { mapped_reads }

   }

   else if ( params.mapper == "bowtie2" ) {

      pre_mapper
         | bowtie2_align
         | set { mapped_reads }

   } 
   
   else if ( params.mapper == "minimap2" ) {

      minimap_align(
         pre_mapper,
         Channel.value( params.nanopore ),

      )
         | set { mapped_reads }

   }

   mapped_reads.main
      | remove_multimappers
      | combine( fetch_UMIcollapse.out )  // sample_id, [bam_bai], umicollapse_repo
      | UMIcollapse  // sample_id, dedup_bam

   UMIcollapse.out.main
      | (
         SAMtools_stats
         & SAMtools_coverage
         & SAMtools_flagstat
      )
   SAMtools_stats.out | plot_bamstats

   UMIcollapse.out.main
      .combine( genome_gff, by: 0 )  // sample_id, dedup_bam, gff
      .unique()
      .set { collapsed }

   featurecounts(
      collapsed,
      Channel.value( !params.nanopore ),
      Channel.value( params.strand ),
      Channel.value( params.ann_type ),
      Channel.value( params.label ),
   )

   collapsed
      .map { tuple( it[0], it[2] ) }
      .unique()
      | gff2bed

   collapsed
      .map { tuple( it[0], it[1] ) }
      .unique()
      .combine( 
         gff2bed.out,
         by: 0,
      )
      | gene_body_coverage
   
   featurecounts.out.main
   | bam2table
   | (
      counts_per_gene
      & counts_per_cell
      & counts_per_gene_per_cell
   )
   
   // UMItools_count(
   //    featurecounts.out.main,
   //    Channel.value( !params.nanopore ),
   // ) // sample_id, counts

   counts_per_gene_per_cell.out
      .combine( featurecounts.out.table, by: 0 )  // sample_id, umitools_counts, featurecounts_counts
      | join_featurecounts_UMItools
   join_featurecounts_UMItools.out
      | ( 
         plot_UMI_distributions 
         & count_genomes_per_cell 
         & plot_counts_per_gene
         & plot_counts_per_cell
         & build_AnnData
      )

   build_AnnData.out.main | filter_AnnData
   filter_AnnData.out.main | cluster_cells

   trim_logs
      .map { it[1] }
      .concat(
         fastQC.out.multiqc_logs,
         mapped_reads.logs,
         featurecounts.out.logs,
         gene_body_coverage.out.logs,
         SAMtools_coverage.out.map { it[1] },
         SAMtools_flagstat.out.map { it[1] },
         SAMtools_stats.out.map { it[1] },
      )
      .flatten()
      .unique()
      .collect()
      | multiQC

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