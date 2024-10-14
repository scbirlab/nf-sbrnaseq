# Illumina paired-end sbRNA-seq pipeline

Nextflow pipeline to process *already demultiplexed* Illumina paired-end FASTQ files from multiple bacterial samples into a gene $\times$ cell count table.

## Processing steps

For each set of barcodes:

1. Generate all combinations to make a cell barcode whitelist.
2. For each cell barcode, generate all 1-error variants for parsing by `umi-tools`.

For each reference genome:

1. Download the genome FASTA and GFF.

For each sample:

1. Trim reads to adapters using `cutadapt`. Reads from sbRNA-seq will be flanked by sequences containing cell barcodes and UMIs, so this step trims any extra sequences either side of these flanking sequences.
2. Extract cell barcodes and UMIs using `umi-tools`.
3. Map to genome FASTA using `bowtie2`.
4. Deduplicate mapped reads using `umitools`.
5. Count deduplicated reads per gene using `featureCounts`.
6. Count deduplicated reads per gene per cell using `umi-tools`.

### Downstream analyses [work in progress, not yet implemented]

1. De novo transcriptome assembly with Trinity.
2. Isoform analysis with RSEM.
3. Plotting and visualisation.

### Other steps

1. Get FASTQ quality metrics with `fastqc`.
2. Compile the logs of processing steps into an HTML report with `multiqc`.

## Requirements

### Software

You need to have Nextflow and `conda` installed on your system.

### First time using Nextflow?

If it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet)  and, optionally, a [`nextflow.config` file](#inputs) in the directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-sbrnaseq
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which will not be automatically updated.  If you want to ensure that you're running a version of 
the pipeline, use the `-r <version>` flag. For example,

```bash 
nextflow run scbirlab/nf-sbrnaseq -r v0.0.1
```

A list of versions can be found by running `nextflow info scbirlab/nf-sbrnaseq`.

For help, use `nextflow run scbirlab/nf-sbrnaseq --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV containing sample IDs matched with FASTQ filenames, genome information, and adapter sequences.
- `fastq_dir`: path to directory containing the FASTQ files

The following parameters have default values can be overridden if necessary.

- `inputs = "inputs"`: path to directory containing files referenced in the `sample_sheet`, such as lists of guide RNAs.
- `output = "outputs"`: path to directory to put output files 
- `trim_qual = 5` : For `cutadapt`, the minimum Phred score for trimming 3' calls
- `min_length = "9:38"` : For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
- `strand = 1` : For `featureCounts`, the strandedness of RNA-seq. `1` for forward, `2` for reverse.
- `ann_type = 'gene'` : For `featureCounts`, features from GFF column 3 to use for counting
- `label = 'Name'` : For `featureCounts`, one or more (comma-separated) fields from column 9 of GFF for labeling counts

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
   
    sample_sheet = "/path/to/sample_sheet.csv"
    fastq_dir = "/path/to/fastqs"

    trim_qual = 15
}
```

Alternatively, you can provide these on the command line:

```bash
nextflow run scbirlab/nf-sbrnaseq \
    --sample_sheet /path/to/sample_sheet.csv \
    --fastq_dir /path/to/fastqs \
    --trim_qual 15
``` 

### Sample sheet

The sample sheet is a CSV file providing information about each sample: which FASTQ files belong 
to it, the reference genome accession number, adapters to be trimmed off, and the UMI and cell barcode scheme for each sample.

The file must have a header with the column names below, and one line per sample to be processed.

- `sample_id`: the unique name of the sample
- `genome_id`: The [NCBI assembly accession](https://www.ncbi.nlm.nih.gov/datasets/genome/) number for the organism that the guide RNAs are targeting. This number starts with "GCF_" or "GCA_".
- `fastq_pattern`: the search glob to find FASTQ files for each sample in `fastq_dir`. The pipleine will look for files matching `<fastq_dir>/*<fastq_pattern>*`, and should match only two files, corresponding to paired reads.
- `adapter_read1_3prime`: the 3' adapter on the forward read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). The adapter itself and sequences downstream will be removed.
- `adapter_read2_3prime`:  the 3' adapter on the reverse read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). The adapter itself and sequences downstream will be removed.
- `adapter_read1_5prime`: the 5' adapter on the forward read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequences _upstream_ will be removed, but the adapters themselves will be retained.
- `adapter_read2_5prime`:  the 5' adapter on the reverse read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequences _upstream_ will be removed, but the adapters themselves will be retained.
- `umi_read1`: the cell barcode and UMI pattern in [`umitools` regex format](https://umi-tools.readthedocs.io/en/latest/regex.html#regex-regular-expression-mode) for the forward read
- `umi_read2`: the cell barcode and UMI pattern in [`umitools` regex format](https://umi-tools.readthedocs.io/en/latest/regex.html#regex-regular-expression-mode) for the reverse read

Here is an example of the sample sheet:

| sample_id | fastq_pattern | genome_id | adapter_read1_3prime | adapter_read2_3prime | adapter_read1_5prime | adapter_read2_5prime | umi_read1 | umi_read2 |
| --------- | --------- | ------------- | ------------- | --------- | --------- | --------- | --------- | ------------- |
| EcoHX1 | G5512A22_R | GCF_904425475.1 | CAGN{6}G{3} | N{7}N{8}TTATTATA | TATAATAAN{8}N{7} | C{3}N{6}CTG | ^(?P<discard_1>.{3})(?P<cell_1>.{6}).* | ^(?P<discard_2>.{8})(?P<cell_2>.{8})(?P<umi_1>.{7}).* |
| EcoHX2 | G5512A23_R | GCF_904425475.1 | CAGN{6}G{3} | N{7}N{8}TTATTATA | TATAATAAN{8}N{7} | C{3}N{6}CTG | ^(?P<discard_1>.{3})(?P<cell_1>.{6}).* | ^(?P<discard_2>.{8})(?P<cell_2>.{8})(?P<umi_1>.{7}).* |

## Outputs

Outputs are saved in the same directory as `sample_sheet`. They are organised under three directories:

- `processed`: FASTQ files and logs resulting from trimming and UMI extraction
- `counts`: tables and BAM files corresponding to cell $\times$ gene counts
- `multiqc`: HTML report on processing steps

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-sbrnaseq/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [umitools](https://umi-tools.readthedocs.io/en/latest/index.html)
- [featureCounts](https://subread.sourceforge.net/featureCounts.html)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools](http://www.htslib.org/doc/samtools.html)