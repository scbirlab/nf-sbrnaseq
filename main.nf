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

   TRIM_CUTADAPT.out.main  // sample_id, [reads]
      .combine( umi_ch, by: 0 )  // sample_id, [reads], umis
      .combine( whitelists, by: 0 )  // sample_id, [reads], umis, whitelist
      | UMITOOLS_EXTRACT  // sample_id, [reads]

   UMITOOLS_EXTRACT.out.main
      .combine( genome_idx, by: 0 )  // sample_id, [reads], [genome_idx], genome_acc
      | BOWTIE2_ALIGN  // sample_id, [bam_bai]

   BOWTIE2_ALIGN.out.main
      .combine( DOWNLOAD_UMICOLLAPSE.out )  // sample_id, [bam_bai], umicollapse_repo
      | UMICOLLAPSE  // sample_id, dedup_bam

   UMICOLLAPSE.out.main
      .combine( genome_gff, by: 0 )  // sample_id, dedup_bam, gff
      | FEATURECOUNTS  // sample_id, [count_bam_bai]

   FEATURECOUNTS.out.main
      | UMITOOLS_COUNT   // sample_id, counts
   UMITOOLS_COUNT.out.main
      | MAKE_COUNT_MATRIX
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
      -d ${params.min_length} -D 10000  \
      -t ${params.ann_type} -g ${params.label} \
      -a ${gff} \
      -R BAM \
      -T ${task.cpus} \
      --countReadPairs \
      --verbose \
      --extraAttributes gene_biotype,locus_tag \
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
      --stdout ${sample_id}.umitools_count.tsv
   rm ${bamfile.getBaseName()}.sorted.bam
   """
}


process MAKE_COUNT_MATRIX {
   tag "${sample_id}"
   label "big_mem"

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
   label "big_mem"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( counts_table )

   output:
   tuple val( sample_id ), path( "*.h5ad" ), path( "*.png" )

   script:
   """
   #!/usr/bin/env python

   import anndata as ad
   from carabiner import print_err
   from carabiner.mpl import figsaver, grid
   import pandas as pd
   import numpy as np
   import scanpy as sc

   from scipy.sparse import csr_matrix
   from pandas.api.types import CategoricalDtype

   CELL_ID = "cell_barcode"
   GENE_ID = "gene_id"
   MIN_GENES_PER_CELL = 20
   MIN_CELLS_PER_GENE = 1

   print_err("Loading ${counts_table} as a sparse matrix...")
   df = pd.read_csv("${counts_table}", sep="\\t", low_memory=False)
   cell_cat = CategoricalDtype(sorted(df[CELL_ID].unique()), ordered=True)
   gene_cat = CategoricalDtype(sorted(df[GENE_ID].unique()), ordered=True)

   cell_idx = df[CELL_ID].astype(cell_cat).cat.codes
   gene_idx = df[GENE_ID].astype(gene_cat).cat.codes
   matrix = csr_matrix(
      (df["umi_count"].astype(float), (cell_idx, gene_idx)),
      shape=(cell_cat.categories.size, gene_cat.categories.size)
   )
   matrix_df = pd.DataFrame.sparse.from_spmatrix(
      matrix,
      index=cell_cat.categories,
      columns=gene_cat.categories,
   )
   adata = ad.AnnData(matrix_df)
   print_err(adata)
   adata.uns["sample_id"] = "${sample_id}"

   print_err("Annotating gene features...")
   gene_ann_columns = ["locus_tag", "gene_biotype", "featurecounts_count", "Length"]
   gene_ann_df = df[[GENE_ID] + gene_ann_columns].drop_duplicates().set_index(GENE_ID)
   gene_ann_df = gene_ann_df.loc[adata.var_names,:]
   for col in gene_ann_columns:
      adata.var[col.casefold()] = gene_ann_df[col]

   adata.var["ribo"] = (adata.var["gene_biotype"] == "rRNA")
   adata.var["trna"] = (adata.var["gene_biotype"] == "tRNA")
   adata.var["ncrna"] = (adata.var["gene_biotype"] == "ncRNA")
   print_err(adata)

   print_err("Calculating QC metrics...")
   sc.pp.calculate_qc_metrics(
      adata, 
      qc_vars=["ribo", "trna", "ncrna"], 
      inplace=True, 
      percent_top=[20], 
      log1p=True,
   )
   print_err(adata)
   print_err(f"Filtering cells with {MIN_GENES_PER_CELL=}...")
   sc.pp.filter_cells(
      adata, 
      min_genes=MIN_GENES_PER_CELL,
   )
   print_err(adata)
   print_err(f"Filtering genes with {MIN_CELLS_PER_GENE=}...")
   sc.pp.filter_genes(
      adata, 
      min_cells=MIN_CELLS_PER_GENE,
   )
   print_err(adata)
   print_err("Applying PCA...")
   sc.pp.pca(
      adata, 
   )
   print_err("Getting nearest neighbours...")
   sc.pp.neighbors(
      adata, 
      n_neighbors=15, 
      metric='cosine',
   )
   print_err("Carrying out Leiden clustering...")
   sc.tl.leiden(
      adata, 
      key_added="leiden_res1", 
      resolution=1.,
      flavor='igraph',
      directed=False,
   )
   print_err(f"Found {adata.obs['leiden_res1'].cat.categories.size=} Leiden clusters!")
   print_err("Embedding as UMAP...")
   sc.tl.umap(
      adata, 
   )
   print_err("Saving as ${sample_id}.h5ad...")
   adata.write("${sample_id}.h5ad", compression="gzip")

   colors_to_plot = [
      "total_counts", 
      "pct_counts_ribo", 
      "pct_counts_trna",
      "pct_counts_ncrna",
      "leiden_res1",
   ]
   print_err("Making UMAP plots...")
   fig, axes = grid(ncol=len(colors_to_plot), aspect_ratio=1.1)
   for ax, col in zip(axes, colors_to_plot):
      colors = adata.obs[col]
      if isinstance(colors, pd.Series) and colors.dtype == "category":
         colors = [f"C{i}" for i in colors.cat.codes]
      sc = ax.scatter(
         adata.obsm['X_umap'][:,0],
         adata.obsm['X_umap'][:,1],
         c=colors,
         cmap="cividis",
         s=.5,
         alpha=.7,
         norm="log" if col == "total_counts" else None,
      )
      if isinstance(colors, pd.Series) and colors.dtype != "category":
         cb = fig.colorbar(sc, ax=ax)
      ax.set(title=col, xlabel="UMAP_1", ylabel="UMAP_2")
   print_err("Saving UMAP plots as ${sample_id}.umap...")
   figsaver()(
      fig=fig,
      name='${sample_id}.umap',
   )
   print_err("Done!")
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
   tuple val( sample_id ), path( "*.umis-and-cells-per-gene.tsv" ), emit: cells_per_gene
   tuple val( sample_id ), path( "*.genes-per-cell.tsv" ), path( "*.genes-per-biotype-per-cell.tsv" ), emit: genes_per_cell

   script:
   """
   set -x
   grep -v '^#' ${featurecounts_table} | tr \$'\\t' , > fc0.csv
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
   figsaver()(
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
   figsaver()(
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
      figsaver()(
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