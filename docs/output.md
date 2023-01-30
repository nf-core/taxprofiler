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
- [Bracken](#bracken) - Taxonomic classifier using k-mers and abundance estimations
- [Kraken2](#kraken2) - Taxonomic classifier using exact k-mer matches
- [KrakenUniq](#krakenuniq) - Taxonomic classifier that combines the k-mer-based classification and the number of unique k-mers found in each species
- [Centrifuge](#centrifuge) - Taxonomic classifier that uses a novel indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index.
- [Kaiju](#kaiju) - Taxonomic classifier that finds maximum (in-)exact matches on the protein-level.
- [Diamond](#diamond) - Sequence aligner for protein and translated DNA searches.
- [MALT](#malt) - Sequence alignment and analysis tool designed for processing high-throughput sequencing data, especially in the context of metagenomics
- [MetaPhlAn3](#metaphlan3) - Genome-level marker gene based taxonomic classifier
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

[fastp](https://github.com/OpenGene/fastp) is a FASTQ pre-processing tool for quality control, trimmming of adapters, quality filtering and other features.

It is used in nf-core/taxprofiler for adapter trimming of short-reads.

<details markdown="1">
<summary>Output files</summary>

- `fastp`
  - `<sample_id>.fastp.fastq.gz`: File with the trimmed unmerged fastq reads.
  - `<sample_id>.merged.fastq.gz`: File with the reads that were successfully merged.

</details>

By default nf-core/taxprofiler will only provide the `<sample_id>.fastp.fastq.gz` file if fastp is selected. The file `<sample_id>.merged.fastq.gz` will be available in the output folder if you provide the argument ` --shortread_qc_mergepairs` (optionally retaining un-merged pairs when in combination with `--shortread_qc_includeunmerged`). You can change the default value for low complexity filtering by using the argument `--shortread_complexityfilter_fastp_threshold`.

### AdapterRemoval

[AdapterRemoval](https://adapterremoval.readthedocs.io/en/stable/) searches for and removes remnant adapter sequences from High-Throughput Sequencing (HTS) data and (optionally) trims low quality bases from the 3' end of reads following adapter removal. It is popular in the field of palaeogenomics. The output logs are stored in the results folder, and as a part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

- `adapterremoval/`
  - `<sample_id>.settings`: AdapterRemoval log file containing general adapter removal, read trimming and merging statistics
  - `<sample_id>.collapsed.fastq.gz` - read-pairs that merged and did not undergo trimming (only when `--shortread_qc_mergepairs` supplied)
  - `<sample_id>.collapsed.truncated.fastq.gz` - read-pairs that merged underwent quality trimming (only when `--shortread_qc_mergepairs` supplied)
  - `<sample_id>.pair1.truncated.fastq.gz` - read 1 of pairs that underwent quality trimming
  - `<sample_id>.pair2.truncated.fastq.gz` - read 2 of pairs that underwent quality trimming (and could not merge if `--shortread_qc_mergepairs` supplied)
  - `<sample_id>.singleton.truncated.fastq.gz` - orphaned read pairs where one of the pair was discarded
  - `<sample_id>.discard.fastq.gz` - reads that were discarded due to length or quality filtering

</details>

By default nf-core/taxprofiler will only provide the `.settings` file if AdapterRemoval is selected. You will only find the `.fastq` files in the results directory if you provide ` --save_preprocessed_reads` . If this is selected, you may receive different combinations of `.fastq` files for each sample depending on the input types - e.g. whether you have merged or not, or if you're supplying both single- and paired-end reads.

Note that the `.fastq` files may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as complexity filtering, host removal, run merging etc..

### Porechop

[Porechop](https://github.com/rrwick/Porechop) is a tool for finding and removing adapters from Oxford Nanopore reads. Adapters on the ends of reads are trimmed and if a read has an adapter in its middle, it is considered a chimeric and it chopped into separate reads.

<details markdown="1">
<summary>Output files</summary>

- `porechop`
  - `<sample_id>.fastq.gz`: Adapter-trimmed file

</details>

The output logs are saved in the output folder and are part of MultiQC report.You do not normally need to check these manually.
We do **not** recommend using Porechop if you are already trimming the adapters with ONT's basecaller Guppy.

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

Note that the `.fastq` file(s) may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as host removal, run merging etc..

### PRINSEQ++

[PRINSEQ++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) is a C++ implementation of the [prinseq-lite.pl](https://prinseq.sourceforge.net/) program. It can be used to filter, reformat or trim genomic and metagenomic sequence data.

It is used in nf-core/taxprofiler for complexity filtering using different algorithms. This means that it will remove reads with low sequence diversity (e.g. mono- or dinucleotide repeats).

<details markdown="1">
<summary>Output files</summary>

- `prinseqplusplus/`
  - `<sample_id>.log`: log file containing number of reads. Row IDs correspond to: `min_len, max_len, min_gc, max_gc, min_qual_score, min_qual_mean, ns_max_n, noiupac, derep, lc_entropy, lc_dust, trim_tail_left, trim_tail_right, trim_qual_left, trim_qual_right, trim_left, trim_right`
  - `<sample_id>_good_out.fastq.gz`: resulting FASTQ file without low-complexity reads

</details>

By default nf-core/taxprofiler will only provide the `.log` file if PRINSEQ++ is selected as the complexity filtering tool. You will only find the complexity filtered `.fastq` files in your results directory if you supply ` --save_complexityfiltered_reads` .

Note that the `.fastq` file(s) may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as host removal, run merging etc..

### Filtlong

[Filtlong](https://github.com/rrwick/Filtlong) is a quality filtering tool for long reads. It can take a set of small reads and produce a smaller, better subset.

<details markdown="1">
<summary>Output files</summary>

- `filtlong`
  - `<sample_id>_filtered.fastq.gz`: Quality or short read data filtered file
  - `<sample_id>_filtered.log`: log file containing summary statistics

</details>

You can use the filtered `.fastq` for other downstream analyses to reduce repeated preprocessing of files.

We do **not** recommend using Filtlong if you are performing filtering of low quality reads with ONT's basecaller Guppy.

### Bowtie2

[Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes.

It is used with nf-core/taxprofiler to allow removal of 'host' (e.g. human) or other possible contaminant reads (e.g. Phi X) from short-read `.fastq` files prior to profiling.

<details markdown="1">
<summary>Output files</summary>

- `bowtie2/`
  - `<sample_id>.bam`: reads that aligned against the user-supplied reference genome
  - `<sample_id>.bowtie2.log`: log file about the mapped reads
  - `<sample_id>.unmapped.fastq.gz`: the off-target reads from the mapping that is used in downstream steps.

</details>

By default nf-core/taxprofiler will only provide the `.log` file if host removal is turned on. You will only see the mapped (host) reads `.bam` file or the off-target reads in `.fastq` format in your results directory if you provide `--save_hostremoval_mapped` and ` --save_hostremoval_unmapped` respectively.

Note that the `.fastq` file(s) may _not_ always be the 'final' reads that go into taxonomic classification, if you also run other steps such as host removal, run merging etc.. Furthermore, while there is a dedicated section in the MultiQC HTML for Bowtie2, these values are not displayed by default in the General Stats table. Rather, alignment statistics to host genome is reported via samtools stats module in MultiQC report for direct comparison with minimap2 (see below).

### minimap2

[minimap2](https://github.com/lh3/minimap2) is an alignment tool suited to mapping long reads to reference sequences.

It is used with nf-core/taxprofiler to allow removal of 'host' (e.g. human) or other possible contaminant reads from long-read `.fastq` files prior to taxonomic classification/profiling.

<details markdown="1">
<summary>Output files</summary>

- `minimap2`
  - `<sample_id>.bam`: Alignment file in bam format

</details>

By default, nf-core taxprofiler will only provide the `.bam` file if host removal for long reads is turned on.
Note: minimap2 is not yet supported as a module in MultiQC and therefore there is no dedicated section in the MultiQC HTML. Rather, alignment statistics to host genome is reported via samtools stats module in MultiQC report.

### Samtools stats

[Samtools stats](http://www.htslib.org/doc/samtools-stats.html) collects statistics from a SAM, `.bam`, or CRAM alignment file and outputs in a text format.

<details markdown="1">
<summary>Output files</summary>

- `samtoolsstats`
  - `<sample_id>.stats`: File containing samtools stats output

</details>

In most cases you do not need to check this file, as it is rendered in the MultiQC run report.

### Bracken

[Bracken](https://ccb.jhu.edu/software/bracken/) (Bayesian Reestimation of Abundance with Kraken) is a highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample. Braken uses the taxonomy labels assigned by Kraken, a highly accurate metagenomics classification algorithm, to estimate the number of reads originating from each species present in a sample.

> 🛈 The first step of using Bracken requires running Kraken2, therefore the initial results before abundance estimation will be found in `<your_results>/kraken2/<your_bracken_db_name>`.

<details markdown="1">
<summary>Output files</summary>

- `bracken/`
  - `bracken_<db_name>_combined_reports.txt`: combined bracken results as output from Bracken's `combine_bracken_outputs.py` script
  - `<db_name>/`
    - `<sample>_<db_name>.tsv`: TSV file containing per-sample summary of Bracken results with abundance information

</details>

The main taxonomic profiling file from Bracken is the `*.tsv` file. This provides the basic results from Kraken2 but with the corrected abundance information.

### Kraken2

[Kraken](https://ccb.jhu.edu/software/kraken2/) is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences. Kraken examines the k-mers within a query sequence and uses the information within those k-mers to query a database. That database maps -mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  - `<db_name>_combined_reports.txt`: A combined profile of all samples aligned to a given database (as generated by `krakentools`)
  - `<db_name>/`
    - `<sample_id>_<db_name>.classified.fastq.gz`: FASTQ file containing all reads that had a hit against a reference in the database for a given sample
    - `<sample_id>_<db_name>.unclassified.fastq.gz`: FASTQ file containing all reads that did not have a hit in the database for a given sample
    - `<sample_id>_<db_name>.report.txt`: A Kraken2 report that summarises the fraction abundance, taxonomic ID, number of Kmers, taxonomic path of all the hits in the Kraken2 run for a given sample
    - `<sample_id>_<db_name>.classifiedreads.txt`: A list of read IDs and the hits each read had against each database for a given sample

</details>

The main taxonomic classification file from Kraken2 is the `_combined_reports.txt` or `*report.txt` file. The former provides you the broadest over view of the taxonomic classification results across all samples against a single databse, where you get two columns for each sample e.g. `2_all` and `2_lvl`, as well as a summarised column summing up across all samples `tot_all` and `tot_lvl`. The latter gives you the most information for a single sample. The report file is also used for the taxpasta step.

You will only receive the `.fastq` and `*classifiedreads.txt` file if you supply `--kraken2_save_reads` and/or `--kraken2_save_readclassification` parameters to the pipeline.

### KrakenUniq

[KrakenUniq](https://github.com/fbreitwieser/krakenuniq) (formerly KrakenHLL) is an extenson to the fast k-mer-based classification [Kraken](https://github.com/DerrickWood/kraken) with an efficient algorithm for additionally assessing the coverage of unique k-mers found in each species in a dataset.

<details markdown="1">
<summary>Output files</summary>

- `krakenuniq/`
  - `<db_name>/`
    - `<sample_id>_<db_name>.classified.fastq.gz`: FASTQ file containing all reads that had a hit against a reference in the database for a given sample
    - `<sample_id>_<db_name>.unclassified.fastq.gz`: FASTQ file containing all reads that did not have a hit in the database for a given sample
    - `<sample_id>_<db_name>.report.txt`: A Kraken2-style report that summarises the fraction abundance, taxonomic ID, number of Kmers, taxonomic path of all the hits, with an additional column for k-mer coverage, that allows for more accurate distinguishing between false-positive/true-postitive hits
    - `<sample_id>_<db_name>.classifiedreads.txt`: A list of read IDs and the hits each read had against each database for a given sample

</details>

The main taxonomic classification file from KrakenUniq is the `*report.txt` file. This is an extension of the Kraken2 report with the additional k-mer coverage information that provides more information about the accuracy of hits.

> ⚠️ The output system of KrakenUniq can result in other `stdout` or `stderr` logging information being saved in the report file, therefore you must check your report files before downstream use!

You will only receive the `.fastq` and `*classifiedreads.txt` file if you supply `--krakenuniq_save_reads` and/or `--krakenuniq_save_readclassification` parameters to the pipeline.

### Centrifuge

[Centrifuge](https://github.com/DaehwanKimLab/centrifuge) is a taxonomic sequence classifier that uses a Burrows-Wheeler transform and Ferragina-Manzina index for storing and mapping sequences.

<details markdown="1">
<summary>Output files</summary>

- `centrifuge`
  - `<sample_id>.centrifuge.mapped.fastq.gz`: Fastq files containing all mapped reads
  - `<sample_id>.centrifuge.report.txt`: A classification report that summarises the taxonomic ID, the taxonomic rank, length of genome sequence, number of classified and uniquely classified reads
  - `<sample_id>.centrifuge.results.txt`: A file that summarises the classification assignment for a read, i.e read ID, sequence ID, score for the classification, score for the next best classification, number of classifications for this read
  - `<sample_id>.centrifuge.txt`: A Kraken2-style report that summarises the fraction abundance, taxonomic ID, number of k-mers, taxonomic path of all the hits in the centrifuge run for a given sample
  - `<sample_id>.centrifuge.unmapped.fastq.gz`: Fastq file containing all unmapped reads

</details>

The main taxonomic classification files from Centrifuge are the `_combined_reports.txt`, `*report.txt`, `*results.txt` and the `*centrifuge.txt`. The latter is used by the taxpasta step. You will receive the `.fastq` files if you supply `--centrifuge_save_reads`.

### Kaiju

[Kaiju](https://github.com/bioinformatics-centre/kaiju) is a taxonomic classifier that finds maximum exact matches on the protein-level using the Burrows–Wheeler transform.

<details markdown="1">
<summary>Output files</summary>

- `kaiju`
  - `<sample_id>.tsv`: A file that summarises the fraction abundance, taxonomic ID, number of reads and taxonomic names
  - `kaiju_<db_name>_combined_reports.txt`: A combined profile of all samples aligned to a given database (as generated by `kaiju2table`)

</details>

The most summary file is the `*combined_reports.txt` file which summarises results across all samples. However if you wish to look at more precise information about each assignment, check the per-sample file. The default taxonomic rank is `species`. You can provide a different one by updating the argument `--kaiju_taxon_rank`.

### DIAMOND

[DIAMOND](https://github.com/bbuchfink/diamond) is a sequence aligner for translated DNA searches or protein sequences against a protein reference database such as NR. It is a replacement for the NCBI BLAST software tools.It has many key features and it is used as taxonomic classifier in nf-core/taxprofiler.

<details markdown="1">
<summary>Output files</summary>

- `diamond`
  - `<sample_id>.log`: A log file containing stdout information
  - `<sample_id>.sam`: A file in SAM format that contains the aligned reads

</details>

By default you will receive a TSV output. Alternatively, you will receive a `*.sam` file if you provide the parameter `--diamond_save_reads` but in this case no taxonomic classification will be available(!), only the aligned reads in sam format. Note that DIAMOND has many output formats, so depending on your [choice](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options) with ` --diamond_output_format` you will receive the taxonomic information in a different format.

### MALT

[MALT](https://software-ab.cs.uni-tuebingen.de/download/malt) is a fast replacement for BLASTX, BLASTP and BLASTN, and provides both local and semi-global alignment capabilities.

<details markdown="1">
<summary>Output files</summary>

- `malt/`
  - `<db_name>/`
  - `<sample_id>.blastn.sam`: sparse SAM file containing alignments of each hit
  - `<sample_id>.megan`: summary file that can be loaded into the [MEGAN6](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/) interactive viewer. Generated by MEGAN6 companion tool `rma2info`
  - `<sample_id>.rma6`: binary file containing all alignments and taxonomic information of hits that can be loaded into the [MEGAN6](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/) interactive viewer
  - `<sample_id>.txt.gz`: text file containing taxonomic IDs and read counts against each taxon. Generated by MEGAN6 companion tool `rma2info`

</details>

The main output of MALT is the `.rma6` file format, which can be only loaded into MEGAN and it's related tools. We provide the `rma2info` text files for improved compatibility with spreadsheet programs and other programmtic data manipulation tools, however this has only limited information compared to the 'binary' RMA6 file format (the `txt` file only contains taxonomic ID and count, whereas RMA6 has taxonomic lineage information).

You will only receive the `.sam` and `.megan` files if you supply `--malt_save_reads` and/or `--malt_generate_megansummary` parameters to the pipeline.

### MetaPhlAn3

[MetaPhlAn3](https://github.com/biobakery/metaphlan) is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level resolution via marker genes.

<details markdown="1">
<summary>Output files</summary>

- `metaphlan3/`
  - `metaphlan3_<db_name>_combined_reports.txt`: A combined profile of all samples aligned to a given database (as generated by `metaphlan_merge_tables`)
  - `<db_name>/`
    - `<sample_id>.biom`: taxonomic profile in BIOM format
    - `<sample_id>.bowtie2out.txt`: BowTie2 alignment information (can be re-used for skipping alignment when re-running MetaPhlAn3 with different parameters)
    - `<sample_id>_profile.txt`: MetaPhlAn3 taxonomic profile including abundance estimates

</details>

The main taxonomic profiling file from MetaPhlAn3 is the `*_profile.txt` file. This provides the abundance estimates from MetaPhlAn3 however does not include raw counts by default.

### mOTUs

[mOTUS](https://github.com/motu-tool/mOTUs) is a taxonomic profiler that maps reads to a unique marker specific database and estimates the relative abundance of known and unknown species.

<details markdown="1">
<summary>Output files</summary>

- `motus`
  - `<sample_id>.log`: A log file that contains summary statistics
  - `<sample_id>.out`: A classification file that summarises taxonomic identifiers, by default at the rank of mOTUs (i.e., species level), and their relative abundances in the profiled sample.
  - `motus_<db_name>_combined_reports.txt`: A combined profile of all samples aligned to a given database (as generated by `motus_merge`)

</details>

Normally `*_combined_reports.txt` is the most useful file for downstream analyses, but the per sample `.out` file can provide additional more specific information. By default, nf-core/taxprofiler is providing a column describing NCBI taxonomic ID as this is used in the taxpasta step. You can disable this column by activating the argument `--motus_remove_ncbi_ids`.
You will receive the relative abundance instead of read counts if you provide the argument `--motus_use_relative_abundance`.

### Krona

[Krona](https://github.com/marbl/Krona) allows the exploration of (metagenomic) hierarchical data with interactive zooming, multi-layered pie charts.

Krona charts will be generated by the pipeline for supported tools (Kraken2, Centrifuge, Kaiju, and MALT)

<details markdown="1">
<summary>Output files</summary>

- `krona/`
  - `<tool_name>_<db_name>.html`: per-tool/per-database interactive HTML file containing hierarchical piecharts

</details>

The resulting HTML files can be loaded into your web browser for exploration. Each file will have a dropdown to allow you to switch between each sample aligned against the given database of the tool.

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

All tools in taxprofiler supported by MultiQC will have a dedicated section showing summary statistics of each tool based on information stored in log files. Note that the 'General Stats' table by default will only show statistics referring to pre-processing steps, and will not display possible values from each classifier/profiler, unless turned on by the user within the 'Configure Columns' menu or via a custom MultiQC config file (`--multiqc_config`)

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
