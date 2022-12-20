# nf-core/taxprofiler: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [falco](#falco) - Alternative to FastQC for raw read QC
- [fastp](#fastp) - Adapter trimming for Illumina data
- [AdapterRemoval](#adapterremoval) - Adapter trimming for Illumina data
- [Porechop](#porechop) - Adapter removal for Oxford Nanopore data
- [BBDuk](#bbduk) - Quality trimming and filtering for Illumina data
- [PRINSEQ++](#prinseq++) - Quality trimming and filtering for Illunina data
- [Filtlong](#filtlong) - Quality trimming and filtering for Nanopore data
- [Bowtie2](#bowtie2) - Host removal for Illumina reads
- [minimap2](#minimap2) - Host removal for Nanopore reads
- [samtoolsstats](#samtoolsstats) - Statistics from host removal
- [Kraken2](#kraken2) - Taxonomic classifier using exact k-mer matches
- [KrakenUniq](#krakenuniq) - Taxonomic classifier that combines the k-mer-based classification and the number of unique k-mers found in each species
- [Bracken](#bracken) -
- [Centrifuge](#centrifuge) - Taxonomic classifier that uses a novel indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index.
- [Kaiju](#kaiju) - Taxonomic classifier that finds maximum (in-)exact matches on the protein-level.
- [Diamond](#diamond) - Sequence aligner for protein and translated DNA searches.
- [mOTUs](#motus) - Tool for marker gene-based OTU (mOTU) profiling.
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### fastp

fastp can automatically detect adapter sequences for Illumina data.

<details markdown="1">
<summary>Output files</summary>

- `fastp`

</details>

### AdapterRemoval

[AdapterRemoval](https://adapterremoval.readthedocs.io/en/stable/) searches for and removes remnant adapter sequences from High-Throughput Sequencing (HTS) data and (optionally) trims low quality bases from the 3' end of reads following adapter removal. It is popular in the field of palaeogenomics. The output logs are stored in the results folder, and as a part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

- `adapterremoval/`
  - `<sample_id>.settings`: AdapterRemoval log file containing general adapter removal, read trimming and merging statistics
  - `<sample_id>.collapsed.fastq.gz` - read-pairs that merged and did not undergo trimming (only when  `--shortread_qc_mergepairs` supplied)
  - `<sample_id>.collapsed.truncated.fastq.gz`  - read-pairs that merged underwent quality trimming  (only when  `--shortread_qc_mergepairs` supplied)
  - `<sample_id>.pair1.truncated.fastq.gz`  - read 1 of pairs that underwent quality trimming
  - `<sample_id>.pair2.truncated.fastq.gz`  - read 2 of pairs that underwent quality trimming  (and could not merge if  `--shortread_qc_mergepairs` supplied)
  - `<sample_id>.singleton.truncated.fastq.gz` - orphaned read pairs where one of the pair was discarded
  - `<sample_id>.discard.fastq.gz` - reads that were discarded due to length or quality filtering

</details>

By default nf-core/taxprofiler will only provide the `.settings` file if AdapterRemoval is selected. You will only find the FASTQ files in the results directory if you provide ` --save_preprocessed_reads` . If this is selected, you may recieve different combinations of FASTQ files for each sample depending on the input types - e.g. whether you have merged or not, or if you're supplying both single- and paired-end reads.

Note that the FASTQ files may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as complexity filtering, host removal, run merging etc..

### Porechop

<details markdown="1">
<summary>Output files</summary>

- `porechop`
  - `<sample_id>.fastq.gz`

</details>

### BBDuk

[BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) stands for Decontamination Using Kmers. BBDuk was developed to combine most common data-quality-related trimming, filtering, and masking operations into a single high-performance tool.

It is used in nf-core/taxprofiler for complexity filtering using different algorithms. This means that it will remove reads with low sequence diversity (e.g. mono- or dinucleotide repeats).

<details markdown="1">
<summary>Output files</summary>

- `bbduk/`
  - `<sample_id>.bbduk.log`: log file containing filtering statistics
  - `<sample_id>.fastq.gz`: resulting FASTQ file without low-complexity reads

</details>

By default nf-core/taxprofiler will only provide the `.log` file if BBDuk is selected as the complexity filtering tool. You will only find the complexity filtered reads in your results directory if you provide ` --save_complexityfiltered_reads` .

Note that the FASTQ file(s) may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as host removal, run merging etc..

### PRINSEQ++

[PRINSEQ++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) is a C++ implementation of the [prinseq-lite.pl](https://prinseq.sourceforge.net/) program. It can be used to filter, reformat or trim genomic and metagenomic sequence data.

It is used in nf-core/taxprofiler for complexity filtering using different algorithms. This means that it will remove reads with low sequence diversity (e.g. mono- or dinucleotide repeats).

<details markdown="1">
<summary>Output files</summary>

- `prinseqplusplus/`
  - `<sample_id>.log`: log file containing number of reads. Row IDs correspond to: `min_len, max_len, min_gc, max_gc, min_qual_score, min_qual_mean, ns_max_n, noiupac, derep, lc_entropy, lc_dust, trim_tail_left, trim_tail_right, trim_qual_left, trim_qual_right, trim_left, trim_right`
  - `<sample_id>_good_out.fastq.gz`:  resulting FASTQ file without low-complexity reads

</details>

By default nf-core/taxprofiler will only provide the `.log` file if PRINSEQ++ is selected as the complexity filtering tool. You will only find the complexity filtered FASTQ files in your results directory if you supply ` --save_complexityfiltered_reads` .

Note that the FASTQ file(s) may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as host removal, run merging etc..

### Filtlong

<details markdown="1">
<summary>Output files</summary>

- `filtlong`
  - `<sample_id>_filtered.fastq.g`
  - `<sample_id>_filtered.log`

</details>

### Bowtie2

[Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes.

It is used with nf-core/taxprofiler to allow removal of 'host' (e.g. human) or other possible contaminant reads (e.g. Phi X) from the FASTQ files prior to profiling.

<details markdown="1">
<summary>Output files</summary>

- `bowtie2/`
  - `<sample_id>.bam`: reads that aligned against the user-supplied reference genome
  - `<sample_id>.bowtie2.log`: log file about the mapped reads
  - `<sample_id>.unmapped.fastq.gz`: the off-target reads from the mapping that is used in downstream steps.

</details>

By default nf-core/taxprofiler will only provide the `.log` file if host removal is turned on. You will only see the mapped (host) reads BAM file or the off-target reads in FASTQ format in your results directory if you provide `--save_hostremoval_mapped`  and ` --save_hostremoval_unmapped` respectively.

Note that the FASTQ file(s) may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as host removal, run merging etc..

### minimap2

<details markdown="1">
<summary>Output files</summary>

- `minimap2`
  - `<sample_id>.bam`

</details>

### Samtools stats

<details markdown="1">
<summary>Output files</summary>

- `samtoolsstats`
  - `<sample_id>.stats`

</details>

### Kraken2

<details markdown="1">
<summary>Output files</summary>

- `kraken2`
  - `<sample_id>.classified.fastq.gz`
  - `<sample_id>.unclassified.fastq.gz`
  - `<sample_id>.report.txt`
  - `<sample_id>.classifiedreads.txt`

</details>

### KrakenUniq

<details markdown="1">
<summary>Output files</summary>

- `krakenuniq`

  - `<sample_id>.classified.fastq.gz`
  - `<sample_id>.krakenuniq.classified.txt`
  - `<sample_id>.krakenuniq.report.txt`
  - `<sample_id>.unclassified.fastq.gz`

  ## interleaved?

</details>

### Centrifuge

<details markdown="1">
<summary>Output files</summary>

- `centrifuge`
  - `<sample_id>.centrifuge.mapped.fastq.gz`
  - `<sample_id>.centrifuge.report.txt`
  - `<sample_id>.centrifuge.results.txt`
  - `<sample_id>.centrifuge.unmapped.fastq.gz`

</details>

### Kaiju

<details markdown="1">
<summary>Output files</summary>

- `kaiju`
  - `<sample_id>.tsv`

</details>

### Diamond

<details markdown="1">
<summary>Output files</summary>

- `diamond`
  - `<sample_id>.log`
  - `<sample_id>.sam`

</details>

### mOTUs

<details markdown="1">
<summary>Output files</summary>

- `motus`
  - `<sample_id>.log`
  - `<sample_id>.out`

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
