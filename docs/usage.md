# nf-core/taxprofiler: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/taxprofiler/usage](https://nf-co.re/taxprofiler/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/taxprofiler is a pipeline for highly-parallelised taxonomic classification and profiling of shotgun metagenomic data across multiple tools simultaneously. In addition to multiple classification and profiling tools, at the same time it allows you to performing taxonomic classification and profiling across multiple databases and settings per tool, as well as produces standardised output tables to allow immediate cross comparison of results between tools.

In addition to this page, you can find additional usage information on the following pages:

- [Tutorials](usage/tutorials.md)
- [FAQ and Troubleshooting](usage/faq-troubleshooting.md)

## General Usage

To run nf-core/taxprofiler, at a minimum two you require two inputs:

- a sequencing read samplesheet
- a database samplesheet

Both contain metadata and paths to the data of your input samples and databases.

When running nf-core/taxprofiler, every step and tool is 'opt in'. To run a given classifier or profiler you must make sure to supply both a database in your `<database>.csv` and supply `--run_<profiler>` flag to your command. Omitting either will result in the profiling tool not executing.

nf-core/taxprofiler also includes optional pre-processing (adapter clipping, merge running etc.) or post-processing (visualisation) steps. These are also opt in with a `--perform_<step>` flag. In some cases, the pre- and post-processing steps may also require additional files. Please check the parameters tab of this documentation for more information.

Please see the rest of this page for information about how to prepare input samplesheets and databases and how to run Nextflow pipelines. See the [parameters](https://nf-co.re/taxprofiler/parameters) documentation for more information about specific options the pipeline also offers.

## Samplesheet inputs

nf-core/taxprofiler can accept as input raw or preprocessed single- or paired-end short-read (e.g. Illumina) FASTQ files, long-read FASTQ files (e.g. Oxford Nanopore), or FASTA sequences (available for a subset of profilers).

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below. Furthermore, nf-core/taxprofiler also requires a second comma-separated file of 3 columns with a header row as in the examples below.

This samplesheet is then specified on the command line as follows:

```console
--input '[path to samplesheet file]' --databases '[path to database sheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate different runs FASTQ files of the same sample before performing profiling, when `--perform_runmerging` is supplied. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
2612,run1,ILLUMINA,2612_run1_R1.fq.gz,,
2612,run2,ILLUMINA,2612_run2_R1.fq.gz,,
2612,run3,ILLUMINA,2612_run3_R1.fq.gz,2612_run3_R2.fq.gz,

```

:::warning
Runs of the same sample sequenced on Illumina platforms with a combination of single and paired-end data will **not** be run-wise concatenated, unless pair-merging is specified. In the example above, `run3` will be profiled independently of `run1` and `run2` if pairs are not merged.
:::

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 6 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data, as well as long-read FASTA files may look something like the one below. This is for 6 samples, where `2612` has been sequenced twice.

```console
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
2611,ERR5766174,ILLUMINA,,,/<path>/<to>/fasta/ERX5474930_ERR5766174_1.fa.gz
2612,ERR5766176,ILLUMINA,/<path>/<to>/fastq/ERX5474932_ERR5766176_1.fastq.gz,/<path>/<to>/fastq/ERX5474932_ERR5766176_2.fastq.gz,
2612,ERR5766180,ILLUMINA,/<path>/<to>/fastq/ERX5474936_ERR5766180_1.fastq.gz,,
2613,ERR5766181,ILLUMINA,/<path>/<to>/fastq/ERX5474937_ERR5766181_1.fastq.gz,/<path>/<to>/fastq/ERX5474937_ERR5766181_2.fastq.gz,
ERR3201952,ERR3201952,OXFORD_NANOPORE,/<path>/<to>/fastq/ERR3201952.fastq.gz,,
```

:::warning
Input FASTQ and FASTA files _must_ be gzipped.
:::

:::warning
While one can include both short-read and long-read data in one run, we recommend that you split these across _two_ pipeline runs and database sheets (see below). This will allow classification optimisation for each data type, and make MultiQC run-reports more readable (due to run statistics having vary large number differences).
:::

| Column                | Description                                                                                                                                                                                               |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | Unique sample name [required].                                                                                                                                                                            |
| `run_accession`       | Run ID or name unique for each (pairs of) file(s) .Can also supply sample name again here, if only a single run was generated [required].                                                                 |
| `instrument_platform` | Sequencing platform reads generated on, selected from the EBI ENA [controlled vocabulary](https://www.ebi.ac.uk/ena/portal/api/controlledVocab?field=instrument_platform) [required].                     |
| `fastq_1`             | Path or URL to sequencing reads or for Illumina R1 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if data in FASTA is specified. Cannot be combined with `fasta`. |
| `fastq_2`             | Path or URL to Illumina R2 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if single end data. Cannot be combined with `fasta`.                                    |
| `fasta`               | Path or URL to long-reads or contigs in FASTA format. GZipped compressed files accepted. Can be left empty if data in FASTA is specified. Cannot be combined with `fastq_1` or `fastq_2`.                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Full database sheet

nf-core/taxprofiler supports multiple databases being classified/profiled against in parallel for each tool.

Databases can be supplied either in the form of a compressed `.tar.gz` archive of a directory containing all relevant database files or the path to a directory on the filesystem.

:::warning
nf-core/taxprofiler does not provide any databases by default, nor does it currently generate them for you. This must be performed manually by the user. See bottom of this section for more information of the expected database files, or the [building custom database](usage/tutorials#retrieving-databases-or-building-custom-databases) tutorials.
:::

The pipeline takes the paths and specific classification/profiling parameters of the tool of these databases as input via a four column comma-separated sheet.

:::warning
To allow user freedom, nf-core/taxprofiler does not check for mandatory or the validity of non-file database parameters for correct execution of the tool - excluding options offered via pipeline level parameters! Please validate your database parameters (cross-referencing [parameters](https://nf-co.re/taxprofiler/parameters), and the given tool documentation) before submitting the database sheet! For example, if you don't use the default read length - Bracken will require `-r <read_length>` in the `db_params` column.
:::

An example database sheet can look as follows, where 7 tools are being used, and `malt` and `kraken2` will be used against two databases each.

`kraken2` will be run twice even though only having a single 'dedicated' database because specifying `bracken` implies first running `kraken2` on the `bracken` database, as required by `bracken`.

```console
tool,db_name,db_params,db_path
malt,malt85,-id 85,/<path>/<to>/malt/testdb-malt/
malt,malt95,-id 90,/<path>/<to>/malt/testdb-malt.tar.gz
bracken,db1,;-r 150,/<path>/<to>/bracken/testdb-bracken.tar.gz
kraken2,db2,--quick,/<path>/<to>/kraken2/testdb-kraken2.tar.gz
krakenuniq,db3,,/<path>/<to>/krakenuniq/testdb-krakenuniq.tar.gz
centrifuge,db1,,/<path>/<to>/centrifuge/minigut_cf.tar.gz
metaphlan,db1,,/<path>/<to>/metaphlan/metaphlan_database/
motus,db_mOTU,,/<path>/<to>/motus/motus_database/
ganon,db1,,/<path>/<to>/ganon/test-db-ganon.tar.gz
kmcp,db1,;-I 20,/<path>/<to>/kmcp/test-db-kmcp.tar.gz
```

:::warning
For Bracken and KMCP, which are two step profilers, nf-core/taxprofiler has a special way of passing parameters to each steps!

For Bracken, if you wish to supply any parameters to both the Kraken or Bracken steps or just the Bracken step, you **must** have a _semi-colon_ `;` list in the `db_params` column. This allows you to specify the Kraken2 parameters before and Bracken parameters after the `;`. This is particularly important if you supply a Bracken database with a non-default read length parameter. If you do not have any parameters to specify, you can leave this column empty. If you wish to provide settings to _just_ the Kraken2 step of the Bracken profiling, you can supply a normal string to the column without a semi-colon. If you wish to supply parameters to only Bracken (and keep default Kraken2 parameters), then you supply a string to the column starting with `;` and the Bracken parameters _after_.

Similiarly, for KMCP, if you want to supply parameters for both the first (KMCP search) and the _second step_ (KMCP profile) steps, you **must** have a _semi-colon_ separated`;` list in `db_params`. If you wish to provide parameters to just KMCP search, you do not need the `;`. If you want to supply parameters to just KMCP profile (and keep search parameters at default), then you must start the string with `;` and the KMCP profile parameters come _after_ the semi colon. If you do not wish to modify any parameters, you can leave the column empty (i.e. the `;` is not necessary).

This allows you to specify the KMCP search and the KMCP profile parameters, separated by `;`. If you do not have any parameters to specify, you can leave this as empty.
:::

Column specifications are as follows:

| Column      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tool`      | Taxonomic profiling tool (supported by nf-core/taxprofiler) that the database has been indexed for [required]. Please note that `bracken` also implies running `kraken2` on the same database.                                                                                                                                                                                                                                                                                                                                                                                                              |
| `db_name`   | A unique name per tool for the particular database [required]. Please note that names need to be unique across both `kraken2` and `bracken` as well, even if re-using the same database.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `db_params` | Any parameters of the given taxonomic classifier/profiler that you wish to specify that the taxonomic classifier/profiling tool should use when profiling against this specific database. Can be empty to use taxonomic classifier/profiler defaults. Must not be surrounded by quotes [required]. We generally do not recommend specifying parameters here that turn on/off saving of output files or specifying particular file extensions - this should be already addressed via pipeline parameters. For Bracken databases, must at a minimum contain a `;` separating Kraken2 from Bracken parameters. |
| `db_path`   | Path to the database. Can either be a path to a directory containing the database index files or a `.tar.gz` file which contains the compressed database directory with the same name as the tar archive, minus `.tar.gz` [required].                                                                                                                                                                                                                                                                                                                                                                       |

:::tip
You can also specify the same database directory/file twice (ensuring unique `db_name`s) and specify different parameters for each database to compare the effect of different parameters during classification/profiling.
:::

nf-core/taxprofiler will automatically decompress and extract any compressed archives for you.

The (uncompressed) database paths (`db_path`) for each tool are expected to contain:

- [**Bracken**:](usage/tutorials.md#bracken-custom-database) output of the combined `kraken2-build` and `bracken-build` process.
- [**Centrifuge**:](usage/tutorials.md#centrifuge-custom-database) output of `centrifuge-build`.
- [**DIAMOND**:](usage/tutorials.md#diamond-custom-database) output of `diamond makedb`.
- [**Kaiju**:](usage/tutorials.md#kaiju-custom-database) output of `kaiju-makedb`.
- [**Kraken2**:](usage/tutorials.md#kraken2-custom-database) output of `kraken2-build` command(s).
- [**KrakenUniq**:](usage/tutorials.md#krakenuniq-custom-database) output of `krakenuniq-build` command(s).
- [**MALT**](usage/tutorials.md#malt-custom-database) output of `malt-build`.
- [**MetaPhlAn**:](usage/tutorials.md#metaphlan-custom-database) output of with `metaphlan --install` or downloaded from links on the [MetaPhlAn wiki](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4#customizing-the-database).
- [**mOTUs**:](usage/tutorials.md#motus-custom-database) the directory `db_mOTU/` that is downloaded via `motus downloadDB`.
- [**ganon**:](usage/tutorials.md#ganon-custom-database) output of `ganon build` or `ganon build-custom`.
- [**KMCP**:](usage/tutorials.md#kmcp-custom-database) output of `kmcp index`. Note: `kmcp index` uses the output of an upstream `kmcp compute` step.

:::info
Click the links in the list above for short quick-reference tutorials how to generate custom databases for each tool.
:::

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/taxprofiler --input samplesheet.csv --databases databases.csv --outdir <OUTDIR> -profile docker --run_<TOOL1> --run_<TOOL2>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

When running nf-core/taxprofiler, every step and tool is 'opt in'. To run a given classifier/profiler you must make sure to supply both a database in your `<database>.csv` and supply `--run_<profiler>` flag to your command. Omitting either will result in the classification/profiling tool not executing. If you wish to perform pre-processing (adapter clipping, merge running etc.) or post-processing (visualisation) steps, these are also opt in with a `--perform_<step>` flag. In some cases, the pre- and post-processing steps may also require additional files. Please check the parameters tab of this documentation for more information.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/taxprofiler -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Sequencing quality control

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. nf-core taxprofiler offers [`falco`](https://github.com/smithlabcode/falco) as an drop-in replacement, with supposedly better improvement particularly for long reads.

### Preprocessing Steps

nf-core/taxprofiler offers four main preprocessing steps for preprocessing raw sequencing reads:

- [**Read processing**](#read-processing): adapter clipping and pair-merging.
- [**Redundancy estimation**](#redundancy-estimation): short-read metagenome coverage estimation.
- [**Complexity filtering**](#complexity-filtering): removal of low-sequence complexity reads.
- [**Host read-removal**](#host-read-removal): removal of reads aligning to reference genome(s) of a host.
- [**Run merging**](#run-merging): concatenation of multiple FASTQ chunks/sequencing runs/libraries of a sample.

:::info
You can save the 'final' reads used for classification/profiling from any combination of these steps with `--save_analysis_ready_reads`.
:::

#### Read Processing

Raw sequencing read processing in the form of adapter clipping and paired-end read merging can be activated via the `--perform_shortread_qc` or `--perform_longread_qc` flags.

It is highly recommended to run this on raw reads to remove artifacts from sequencing that can cause false positive identification of taxa (e.g. contaminated reference genomes) and/or skews in taxonomic abundance profiles. If you have public data, normally these should have been corrected for, however you should still check that these steps have indeed been already performed.

There are currently two options for short-read preprocessing: [`fastp`](https://github.com/OpenGene/fastp) or [`adapterremoval`](https://github.com/MikkelSchubert/adapterremoval).

For adapter clipping, you can either rely on the tool's default adapter sequences, or supply your own adapters (`--shortread_qc_adapter1` and `--shortread_qc_adapter2`)
By default, paired-end merging is not activated. In this case paired-end 'alignment' against the reference databases is performed where supported, and if not, supported pairs will be independently classified/profiled. If paired-end merging is activated you can also specify whether to include unmerged reads in the reads sent for classification/profiling (`--shortread_qc_mergepairs` and `--shortread_qc_includeunmerged`).
You can also turn off clipping and only perform paired-end merging, if requested. This can be useful when processing data downloaded from the ENA, SRA, or DDBJ (`--shortread_qc_skipadaptertrim`).
Both tools support length filtering of reads and can be tuned with `--shortread_qc_minlength`. Performing length filtering can be useful to remove short (often low sequencing complexity) sequences that result in unspecific classification and therefore slow down runtime during classification/profiling, with minimal gain.

There is currently one option for long-read Oxford Nanopore processing: [`porechop`](https://github.com/rrwick/Porechop).

For both short-read and long-read preprocessing, you can optionally save the resulting processed reads with `--save_preprocessed_reads`.

#### Redundancy Estimation

Metagenome 'coverage' or sequencing complexity estimations of short-read datasets can be activated with `--perform_shortread_redundancyestimation`.

This turns on checking of read redundancy in a sequencing library using [Nonpareil](https://nonpareil.readthedocs.io/en/latest/), to provide an estimation of whether you have sequenced enough to capture all possible genomes present in your metagenomic sample (with the assumption that once you've sequenced enough, you will keep sequencing PCR amplicons rather than unique reads).

This is only suitable for short-read, and in nf-core/taxprofiler specifically, FASTQ files.

Nonpareil is performed on processed reads (i.e. after fastp or AdapterRemoval). This will run on either the first read of each read pair (as recommended by the authors), or on merged reads.

Before using this tool please note the following caveats:

:::warning

- It is not recommended to run this on deep sequencing data, or very large datasets
  - Nonpareil requires uncompressed FASTQ files, and nf-core/taxprofiler will uncompress these in your working directory, potentially with a extremely large hard-drive footprint.
- Your shortest reads _after_ processing should not go below 24bp (see warning below)
- It is not recommended to keep unmerged (`--shortread_qc_includeunmerged`) reads when using the calculation.

:::info
If you get errors regarding the 'kmer' value is not correct, make sure your shortest reads _after_ processing is not less than 24bp.

If this is the case you will need to specify in a custom config

```nextflow
process {
  withName: NONPAREIL_NONPAREIL {
    ext.args = { "-k <NUMBER>" }
    }
}
```

Where `<NUMBER>` should be at least the shortest read in your library
:::

#### Complexity Filtering

Complexity filtering can be activated via the `--perform_shortread_complexityfilter` flag.

Complexity filtering is primarily a run-time optimisation step. It is not necessary for accurate taxonomic classification/profiling, however it can speed up run-time of each tool by removing reads with low-diversity of nucleotides (e.g. with mono-nucleotide - `AAAAAAAA`, or di-nucleotide repeats `GAGAGAGAGAGAGAG`) that have a low-chance of giving an informative taxonomic ID as they can be associated with many different taxa. Removing these reads therefore saves computational time and resources.

There are currently three options for short-read complexity filtering: [`bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), [`prinseq++`](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus), and [`fastp`](https://github.com/OpenGene/fastp#low-complexity-filter).

There is one option for long-read quality filtering: [`Filtlong`](https://github.com/rrwick/Filtlong)

The tools offer different algorithms and parameters for removing low complexity reads and quality filtering. We therefore recommend reviewing the pipeline's [parameter documentation](https://nf-co.re/taxprofiler/parameters) and the documentation of the tools (see links above) to decide on optimal methods and parameters for your dataset.

You can optionally save the FASTQ output of the run merging with the `--save_complexityfiltered_reads`. If running with `fastp`, complexity filtering happens inclusively within the earlier shortread preprocessing step. Therefore there will not be an independent pipeline step for complexity filtering, and no independent FASTQ file (i.e. `--save_complexityfiltered_reads` will be ignored) - your complexity filtered reads will also be in the `fastp/` folder in the same file(s) as the preprocessed read.

:::warning
For nanopore data: we do not recommend performing any read preprocessing or complexity filtering if you are using ONTs Guppy toolkit for basecalling and post-processing.
:::

#### Host-Read Removal

Removal of possible-host reads from FASTQ files prior classification/profiling can be activated with `--perform_shortread_hostremoval` or `--perform_longread_hostremoval`.

Similarly to complexity filtering, host-removal can be useful for runtime optimisation and reduction in misclassified reads. It is not always necessary to report classification of reads from a host when you already know the host of the sample, therefore you can gain a run-time and computational advantage by removing these prior typically resource-heavy classification/profiling with more efficient methods. Furthermore, particularly with human samples, you can reduce the number of false positives during classification/profiling that occur due to host-sequence contamination in reference genomes on public databases.

nf-core/taxprofiler currently offers host-removal via alignment against a reference genome with Bowtie2 for short reads and minimap2 for long reads, and the use of the unaligned reads for downstream classification/profiling.

You can supply your reference genome in FASTA format with `--hostremoval_reference`. You can also optionally supply a directory containing pre-indexed Bowtie2 index files with `--shortread_hostremoval_index` or a minimap2 `.mmi` file for `--longread_hostremoval_index`, however nf-core/taxprofiler will generate these for you if necessary. Pre-supplying the index directory or files can greatly speed up the process, and these can be re-used.

:::tip
If you have multiple taxa or sequences you wish to remove (e.g., the host genome and then also PhiX - common quality-control reagent during sequencing) you can simply concatenate the FASTAs of each taxa or sequences into a single reference file.
:::

#### Run Merging

For samples that may have been sequenced over multiple runs, or for FASTQ files split into multiple chunks, you can activate the ability to merge across all runs or chunks with `--perform_runmerging`.

For more information how to set up your input samplesheet, see [Multiple runs of the same sample](#multiple-runs-of-the-same-sample).

Activating this functionality will concatenate the FASTQ files with the same sample name _after_ the optional preprocessing steps and _before_ classification/profiling. Note that libraries with runs of different pairing types will **not** be merged and this will be indicated on output files with a `_se` or `_pe` suffix to the sample name accordingly.

You can optionally save the FASTQ output of the run merging with the `--save_runmerged_reads`.

#### Classification and Profiling

The following sections provide tips and suggestions for running the different taxonomic classification and profiling tools _within the pipeline_. For advice and/or guidance whether you should run a particular tool on your specific data, please see the documentation of each tool!

An important distinction between the different tools in included in the pipeline is classification versus profiling. Taxonomic _classification_ is concerned with simply detecting the presence of species in a given sample. Taxonomic _profiling_ involves additionally estimating the _abundance_ of each species.

Note that not all taxonomic classification tools (e.g. Kraken, MALT, Kaiju) performs _profiling_, but all taxonomic profilers (e.g. MetaPhlAn, mOTUs, Bracken) must perform some form of _classification_ prior to profiling.

For advice as to which tool to run in your context, please see the documentation of each tool.

:::note
If you would like to change this behaviour, please contact us on the [nf-core slack](https://nf-co.re/join) and we can discuss this.
:::

Not all tools currently have dedicated tips, suggestions and/or recommendations, however we welcome further contributions for existing and additional tools via pull requests to the [nf-core/taxprofiler repository](https://github.com/nf-core/taxprofiler)!

##### Bracken

You must make sure to also activate Kraken2 to run Bracken in the pipeline.

It is unclear whether Bracken is suitable for running long reads, as it makes certain assumptions about read lengths. Furthermore, during testing we found issues where Bracken would fail on the long-read test data.

Therefore currently nf-core/taxprofiler does not run Bracken on data specified as being sequenced with `OXFORD_NANOPORE` in the input samplesheet.

##### Centrifuge

Centrifuge currently does not accept FASTA files as input, therefore no output will be produced for these input files.

##### DIAMOND

DIAMOND only allows output of a single file format at a time, therefore parameters such `--diamond_save_reads` supplied will result in only aligned reads in SAM format will be produced, no taxonomic profiles will be available. Be aware of this when setting up your pipeline runs, depending on your particular use case.

##### Kaiju

Currently, no specific tips or suggestions.

##### Kraken2

Currently, no specific tips or suggestions.

##### KrakenUniq

Currently, no specific tips or suggestions.

##### MALT

MALT does not support paired-end reads alignment (unlike other tools), therefore nf-core/taxprofiler aligns these as independent files if read-merging is skipped. If you skip merging, you can sum or average the results of the counts of the pairs.

Krona can only be run on MALT output if path to Krona taxonomy database supplied to `--krona_taxonomy_directory`. Therefore if you do not supply the a Krona directory, Krona plots will not be produced for MALT.

##### MetaPhlAn

MetaPhlAn4 is compatible with the MetaPhlAn3 database by adding the `--mpa3` into `db_params` of the `database.csv`.

##### mOTUs

mOTUs currently does not accept FASTA files as input, therefore no output will be produced for these input files.

##### ganon

It is unclear whether ganon is suitable for running long reads - during testing we found issues where ganon would fail on the long-read test data.

Therefore currently nf-core/taxprofiler does not run ganon on data specified as being sequenced with `OXFORD_NANOPORE` in the input samplesheet.

##### KMCP

KMCP is only suitable for short-read metagenomic profiling, with much lower sensitivity on long-read datasets. Therefore, nf-core/taxprofiler does not currently run KMCP on data specified as being sequenced with `OXFORD_NANOPORE` in the input samplesheet.

#### Post Processing

##### Visualisation

nf-core/taxprofiler supports generation of Krona interactive pie chart plots for the following compatible tools.

- Kraken2
- Centrifuge
- Kaiju
- MALT

:::warning
MALT KRONA plots cannot be generated automatically, you must also specify a Krona taxonomy directory with `--krona_taxonomy_directory` if you wish to generate these.
:::

##### Multi-Table Generation

The main multiple-sample table from nf-core/taxprofiler is from a dedicated standalone tool originally developed for the pipeline - [Taxpasta](https://taxpasta.readthedocs.io/en/latest/). When providing `--run_profile_standardisation`, every classifier/profiler and database combination will get a standardised and (if present) multi-sample taxon table in the [`taxpasta/`](https://nf-co.re/taxprofiler/output) directory. These tables are structured in the same way, to facilitate comparison between the results of the classifier/profiler. If multiple samples are provided, `taxpasta merge` will be executed, whereas if only a single sample is provided, `taxpasta standardise` will be executed - the file naming scheme will be the same for both.

In addition to per-sample profiles and standardised Taxpasta output, the pipeline also supports generation of 'native' multi-sample taxonomic profiles (i.e., those generated by the taxonomic profiling tools themselves or additional utility scripts provided by the tool authors), when providing `--run_profile_standardisation` to your pipeline.

These are executed on a per-database level. I.e., you will get a multi-sample taxon table for each database you provide for each tool and will be placed in the same directory as the directories containing the per-sample profiles.

The following tools will produce multi-sample taxon tables:

- **Bracken** (via bracken's `combine_bracken_outputs.py` script)
- **Centrifuge** (via KrakenTools' `combine_kreports.py` script)
- **Kaiju** (via Kaiju's `kaiju2table` tool)
- **Kraken2** (via KrakenTools' `combine_kreports.py` script)
- **MetaPhlAn** (via MetaPhlAn's `merge_metaphlan_tables.py` script)
- **mOTUs** (via the `motus merge` command)
- **ganon** (via the `ganon table` command)

Note that the multi-sample tables from the 'native' tools in each folders are [not inter-operable](https://taxpasta.readthedocs.io/en/latest/tutorials/getting-started/) with each other as they can have different formats and can contain additional and different data. In this case we refer you to use the standardised and merged output from Taxpasta, as described above.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/taxprofiler
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/taxprofiler releases page](https://github.com/nf-core/taxprofiler/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
