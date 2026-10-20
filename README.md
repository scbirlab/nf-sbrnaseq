# sbRNA-seq pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-sbrnaseq/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**scbirlab/nf-sbrnaseq** is a Nextflow pipeline for processing *already demultiplexed* Illumina paired-end or Oxford Nanopore FASTQ files from bacterial single-cell RNA-seq experiments. It produces cell × gene count tables and AnnData objects, and supports optional analysis of RNA base substitutions at nominated sites and genomic RNA 5′-end loci. The resulting measurements can be combined in sparse, per-sample MuData (`.h5mu`) files.

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Credit](#credit)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

**scbirlab/nf-sbrnaseq** carries out the following steps.

For each set of barcodes:

1. Generate the combinations of the supplied barcode rounds to make a cell-barcode whitelist.
2. Optionally generate single-error barcode variants for matching with `umi-tools`.

For each reference genome:

1. Download the reference FASTA and GFF from NCBI. Where multiple assemblies are specified together (for example, a barnyard experiment), combine the references for alignment.
2. If modification-site analysis is enabled, resolve the nominated sites against the reference FASTA and GFF (`ResolveSites`).
3. If 5′-end analysis is enabled, pool molecule ends from samples sharing a reference genome and discover strand-specific clusters with `paraclu` (`ClusterMoleculeEnds`). Select non-overlapping loci and annotate their proximity to same-strand genes with `bedtools closest` (`AnnotateEndLoci`). Clusters are **RNA-end loci**, not automatically transcription start sites.

For each sample:

1. Trim adapter sequences with `cutadapt`, using the Illumina or Nanopore path as appropriate.
2. Extract cell barcodes and UMIs with `umi-tools` and align reads with `STAR`, `bowtie2` or `minimap2`.
3. Apply the configured mapping-quality/multimapping filters; collapse UMIs and create a sorted, indexed, deduplicated molecule-level BAM (`DeduplicateMolecules`).
4. Count genes using `featureCounts` and derive per-cell UMI counts, gene-count tables, mapping QC and coverage plots.
5. If enabled, count observed bases at resolved modification sites per cell (`CountModSites`), including reference, alternative and other base counts, callable depth, and mapping-confidence QC. These are **observed reverse-transcription substitutions**, not a direct measurement of modification stoichiometry.
6. If enabled, extract the selected aligned read boundary for each deduplicated molecule (`ExtractMoleculeEnds`) and assign it by coordinate and strand to at most one selected 5′ locus (`AssignEndLoci`). Ends outside selected loci remain unassigned and are included in assignment QC.
7. Construct the gene-expression AnnData and run existing Scanpy filtering, QC, Leiden clustering and UMAP. Build sparse end5 and modification AnnData objects when enabled (`BuildEnd5AnnData`, `BuildModAnnData`), then assemble a per-sample MuData (`BuildMuData`). Cell identifiers are sample-qualified; cells observed in only some modalities are retained with a `present_in_modality` flag.

### Other steps

1. Assess FASTQ quality with `fastqc`.
2. Generate read-mapping statistics and genome coverage with `samtools`, gene-body coverage with `RSeQC`, and additional coverage plots.
3. Compile processing logs and QC into a `multiqc` HTML report.

### Downstream analyses [work in progress]

Biological interpretation of antibiotic-response heterogeneity, differential 5′-locus usage and modification-associated signatures remains downstream of the pipeline. A 5′-end cluster may arise from transcription initiation, RNA processing or degradation; it is not labelled a TSS by the pipeline.

## Requirements

### Software

Install Nextflow and either Conda, Docker or Singularity/Apptainer. The workflow manages its bioinformatics dependencies through `environment.yml` or the configured container. Important tools include `cutadapt`, `umi-tools`, `STAR`, `bowtie2`, `minimap2`, `samtools`, `featureCounts`, `paraclu`, `bedtools`, `scanpy`, `anndata` and `mudata`.

### First time using Nextflow?

##### Crick users

If you're at the Crick, or your cluster already provides Nextflow and Singularity:

```bash
module load Nextflow Singularity
```

##### Others

Install Nextflow, for example using Conda:

```bash
conda install -c bioconda nextflow
```

You may need to configure `NXF_HOME`:

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

For a persistent shell setting, add the export to your shell initialization file.

## Quick start

Prepare a [sample sheet](#sample-sheet) and run the pipeline with your input FASTQ directory:

```bash
nextflow run scbirlab/nf-sbrnaseq \
    --sample_sheet /path/to/sample_sheet.csv \
    --fastq_dir /path/to/fastqs
```

To try the changes on the development branch, specify `-r fix-multimapping`. For reproducible analyses, pin a release or commit instead of relying on a moving branch:

```bash
nextflow run scbirlab/nf-sbrnaseq -r fix-multimapping \
    --sample_sheet /path/to/sample_sheet.csv \
    --fastq_dir /path/to/fastqs \
    --end5_enabled true
```

Enable or replace the default modification-site rules with `--mod_sites /path/to/sites.tsv`; set `--mod_sites null` in a configuration profile to disable the modification branch. See [Inputs](#inputs) for defaults and read-boundary settings.

For pipeline help, use `nextflow run scbirlab/nf-sbrnaseq --help`. A local `-c nextflow.config` can be used for your settings. Cached copies of remote pipelines are not automatically refreshed; use `-r` to select the intended version.

## Inputs

Required parameters:

- `sample_sheet`: CSV containing sample IDs, FASTQ identifiers, genome accessions and adapter/barcode specifications.
- `fastq_dir`: directory of local FASTQ files, required unless `from_sra = true`.

### General and read-processing parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `sample_sheet` | `null` | Sample sheet path. |
| `fastq_dir` | `null` | Local FASTQ directory. |
| `inputs` | `"inputs"` | Directory of referenced input files, including barcode tables. |
| `outputs` | `"outputs"` | Root output directory. |
| `from_sra` | `false` | Fetch FASTQs by SRA Run accession rather than from local files. |
| `nanopore` | `false` | Enable the Nanopore read-processing path. |
| `unpaired` | `false` | Treat reads as unpaired for applicable downstream operations. |
| `reverse` | `""` | Barcode rounds to reverse-complement when building the whitelist. |
| `allow_cell_errors` | `true` | Allow one error in cell-barcode matching. |
| `whitelist_chunk_size` | `500000000` | Maximum whitelist chunk size. |
| `trim_qual` | `5` | Cutadapt quality threshold. |
| `min_length` | `"9:38"` | Minimum trimmed read lengths as configured for the library. |
| `mapper` | `"star"` | Read aligner: `star`, `bowtie2`, or `minimap2`. |
| `min_mapq` | `20` | Minimum alignment MAPQ for the mapping filter and end extraction. |
| `allow_multimap` | `false` | Allow multimapping alignments in the mapping filter. |
| `max_multimaps` | `20` | Maximum multimapping alignments used in supported mapper configurations. |
| `strand` | `1` | FeatureCounts strandedness: `1` forward, `2` reverse; also used to infer transcript strand for RNA ends. |
| `ann_type` | `"gene"` | Feature type(s) to count from GFF. |
| `label` | `"locus_tag"` | GFF attribute(s) used as featureCounts IDs. |

### Modification-site analysis

| Parameter | Default | Description |
| --- | --- | --- |
| `mod_sites` | `${projectDir}/resources/trna-mod-sites.tsv` | TSV describing nominated sites. Set to `null` to skip site resolution, counting and modification AnnData. |
| `sites_join_on` | `"feature"` | Rule for resolving nominated sites to reference features. |
| `mod_min_mapq` | `20` | Minimum mapping quality for molecule-level site counting. |
| `mod_min_baseq` | `25` | Minimum base quality for site counting. |

`ResolveSites` creates reference-resolved site definitions; `CountModSites` writes per-cell site counts and mapping-confidence statistics. Site-level observations must be interpreted in light of coverage and reference-copy ambiguity. Coverage zero is not the same as an observed reference-only call.

### RNA 5′-end analysis

| Parameter | Default | Description |
| --- | --- | --- |
| `end5_enabled` | `false` | Enable RNA-end extraction, pooled clustering, annotation, assignment and end5 AnnData. |
| `end5_read` | `1` | Read to select: `1` = R1, `2` = R2, `0` = unpaired. Use `0` for truly unpaired reads such as Nanopore alignments. |
| `end5_boundary` | `"5"` | Use the selected sequenced read's `5` or `3` boundary. Configure for your library geometry. |
| `end5_include_softclip` | `false` | Whether to retain alignments soft-clipped at the selected boundary. |
| `end5_min_support` | `5` | Minimum support supplied to Paraclu. |
| `end5_max_width` | `200` | Maximum width of a retained RNA-end locus in nucleotides. |

The end boundary is the *selected sequenced read boundary*, not an intrinsic guarantee of the original RNA end. For the laboratory's R1-based protocol, select R1; other library geometries may need R2 or a different boundary. Positions are reported in 1-based genomic coordinates in the TSV files. Selected loci are non-overlapping on a given reference contig and strand. Clustering uses pooled ends from samples grouped by reference accession; assignment remains per sample and cell.

### Configuration examples

A `nextflow.config` override:

```groovy
params {
    sample_sheet = "/path/to/sample_sheet.csv"
    fastq_dir = "/path/to/fastqs"
    mapper = "star"
    end5_enabled = true
    end5_read = 1
    end5_boundary = "5"
}
```

Or supply settings on the command line:

```bash
nextflow run scbirlab/nf-sbrnaseq -r fix-multimapping \
    --sample_sheet /path/to/sample_sheet.csv \
    --fastq_dir /path/to/fastqs \
    --mapper star \
    --end5_enabled true \
    --end5_read 1 \
    --mod_sites /path/to/modification-sites.tsv
```

### Sample sheet

The sample sheet identifies each sample, its reference genome and the appropriate adapters, UMIs and barcodes. The code uses `genome_accession` as the column name for the NCBI assembly accession (such as `GCF_...` or `GCA_...`), with `+`-separated assemblies for a combined-reference barnyard experiment.

For local FASTQ files (`from_sra = false`):

- `sample_id`: unique sample identifier.
- `genome_accession`: reference assembly accession or `+`-separated combined references.
- `fastq_pattern`: filename fragment/glob used within `fastq_dir`.
- `bc1`, `bc2`, `bc3`: barcode list files, relative to `inputs`.
- `adapter_read1_3prime`, `adapter_read2_3prime`: 3′ adapter patterns for Cutadapt.
- `adapter_read1_5prime`, `adapter_read2_5prime`: 5′ adapter patterns for Cutadapt.
- `umi_read1`, `umi_read2`: cell-barcode and UMI patterns in `umi-tools` regex format.

For SRA input, supply `Run` for each accession. The SRA path prepends Illumina indices to reads; account for those sequences in the adapter and UMI patterns.

For Nanopore samples, the code instead expects `adapter_5prime`, `adapter_3prime` and `umi` in place of the Illumina-specific adapter/UMI columns. Barcode-list and genome fields remain relevant.

Example Illumina sample sheet:

| sample_id | fastq_pattern | genome_accession | bc1 | bc2 | bc3 | adapter_read1_3prime | adapter_read2_3prime | adapter_read1_5prime | adapter_read2_5prime | umi_read1 | umi_read2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EcoHX1 | G5512A22_R | GCF_904425475.1 | bc1.csv | bc2.csv | bc3.csv | CAGN{6}G{3} | N{7}N{8}TTATTATA | TATAATAAN{8}N{7} | C{3}N{6}CTG | ^(?P<discard_1>.{3})(?P<cell_1>.{6}).* | ^(?P<discard_2>.{8})(?P<cell_2>.{8})(?P<umi_1>.{7}).* |

Barcode files must have a header and two columns; the second contains the barcode sequence:

| Barcode ID | Sequence |
| --- | --- |
| barcode-01 | ATCGATCC |
| barcode-02 | TAGCTAGG |

## Outputs

Outputs are written beneath `params.outputs` (default `outputs`). Key directories include:

| Directory | Contents |
| --- | --- |
| `processed/` and extraction/trimming outputs | Intermediate FASTQs and processing logs. |
| `molecules/` | Deduplicated, coordinate-sorted, indexed molecule BAMs (`*.molecules.bam` and `.bai`). |
| `counts/` | Per-cell/per-gene UMI counts and joined count tables. |
| `featurecounts/` | Gene-level featureCounts tables and summaries. |
| `sites/` | Reference-resolved nominated modification sites (`resolved-sites.tsv`). |
| `modifications/` | Per-cell resolved-site base counts and mapping QC (`mod-sites.tsv.gz`, sample-qualified by publishing). |
| `ends/` | Molecule 5′-end tables; Paraclu candidates; selected loci (`loci.tsv`, `loci.bed`); gene-proximity annotations (`annotated.tsv`); cell × locus counts (`counts.tsv`) and assignment QC (`assignment-qc.tsv`), published with sample/reference identifiers. |
| `scanpy/` | Gene-expression `.h5ad`, QC tables and plots. |
| `mudata/` | Sparse `end5.h5ad`, `mod.h5ad` and per-sample `.h5mu` outputs when available. |
| `multiqc/` | Combined HTML QC report. |

### MuData layout

Each sample's `.h5mu` contains an `rna` modality and, if enabled and produced, `end5` and `mod` modalities. All modalities are aligned to the union of their observed cells and use sample-qualified observation identifiers. `present_in_modality` distinguishes cells originating in a modality from zero rows introduced during union reindexing.

- `rna.X`: gene UMI counts (the existing gene AnnData pathway).
- `end5.X`: counts of molecules assigned to each selected genomic RNA-end locus. End-locus genomic coordinates and support are stored in `.var`; extraction/assignment QC is stored in `.obs`.
- `mod.X`: counts of alternative bases at each nominated site. Sparse `.layers` include `ref_count`, `alt_count`, `other_count`, `callable_depth`, `ref_alt_depth` and `count_A`/`count_C`/`count_G`/`count_T`. A zero in `callable_depth` signifies no callable coverage, **not** an observed reference-only call.

The `.h5ad` matrices are not combined into a single expression matrix: they represent different molecular measurements. The existing gene-only Scanpy filtering and clustering remains unchanged.

## Credit

Developed by the SCBIR Lab. Please cite the software and its dependencies when using it in publications; release-specific citation information will be provided when available.

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-sbrnaseq/issues).

## Further help

- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiQC](https://multiqc.info/)
- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [UMItools](https://umi-tools.readthedocs.io/en/latest/index.html)
- [featureCounts](https://subread.sourceforge.net/featureCounts.html)
- [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [minimap2](https://lh3.github.io/minimap2/minimap2.html)
- [samtools](http://www.htslib.org/doc/samtools.html)
- [BEDOPS](https://bedops.readthedocs.io/en/latest/index.html)
- [BEDTools](https://bedtools.readthedocs.io/)
- [Paraclu](https://fantom.gsc.riken.jp/5/suppl/Paraclu/)
- [RSeQC](https://rseqc.sourceforge.net/#)
- [Scanpy](https://scanpy.readthedocs.io/)
- [AnnData](https://anndata.readthedocs.io/)
- [MuData](https://mudata.readthedocs.io/)
