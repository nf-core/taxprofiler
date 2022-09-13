# nf-core/taxprofiler: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/taxprofiler/usage](https://nf-co.re/taxprofiler/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet inputs

nf-core/taxprofiler can accept as input raw or preprocessed single- or paired-end short-read (e.g. Illumina) FASTQ files, long-read FASTQ files (e.g. Oxford Nanopore), or FASTA sequences (available for a subset of profilers).

> ⚠️ Input FASTQ files _must_ be gzipped, while FASTA files may optionally be uncompressed (although this is not recommended)

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below. Furthermother, nf-core/taxprofiler also requires a second comma-separated file of 3 columns with a header row as in the examples below.

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

> ⚠️ Runs of the same sample sequenced on Illumina platforms with a combination of single and paired-end data will **not** be run-wise concatenated, unless pair-merging is specified. In the example above, `run3` will be profiled independently of `run1` and `run2` if pairs are not merged.

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 6 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data, as well as long-read FASTA files may look something like the one below. This is for 6 samples, where `2612` has been sequenced twice.

```console
2611,ERR5766174,ILLUMINA,,,/<path>/<to>/fasta/ERX5474930_ERR5766174_1.fa.gz
2612,ERR5766176,ILLUMINA,/<path>/<to>/fastq/ERX5474932_ERR5766176_1.fastq.gz,/<path>/<to>/fastq/ERX5474932_ERR5766176_2.fastq.gz,
2612,ERR5766180,ILLUMINA,/<path>/<to>/fastq/ERX5474936_ERR5766180_1.fastq.gz,,
2613,ERR5766181,ILLUMINA,/<path>/<to>/fastq/ERX5474937_ERR5766181_1.fastq.gz,/<path>/<to>/fastq/ERX5474937_ERR5766181_2.fastq.gz,
ERR3201952,ERR3201952,OXFORD_NANOPORE,/<path>/<to>/fastq/ERR3201952.fastq.gz,,
```

| Column                | Description                                                                                                                                                                                              |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | Unique sample name [required].                                                                                                                                                                           |
| `run_accession`       | Run ID or name unique for each (pairs of) file(s) .Can also supply sample name again here, if only a single run was generated [required].                                                                |
| `instrument_platform` | Sequencing platform reads generated on, selected from the EBI ENA [controlled vocabulary](https://www.ebi.ac.uk/ena/portal/api/controlledVocab?field=instrument_platform) [required].                    |
| `fastq_1`             | Path or URL to sequencing reads or for Illumina R1 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if data in FASTA is specifed. Cannot be combined with `fasta`. |
| `fastq_2`             | Path or URL to Illumina R2 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if single end data. Cannot be combined with `fasta`.                                   |
| `fasta`               | Path or URL to long-reads or contigs in FASTA format. GZipped compressed files accepted. Can be left empty if data in FASTA is specifed. Cannot be combined with `fastq_1` or `fastq_2`.                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Full database sheet

nf-core/taxprofiler supports multiple databases being profiled in parallel for each tool.
Databases can be supplied either in the form of a compressed `.tar.gz` archive of a directory containing all relevant database files or the path to a directory on the filesystem.
The pipeline takes the locations and specific profiling parameters of the tool of these databases as input via a four column comma-separated sheet.

> ⚠️ nf-core/taxprofiler does not provide any databases by default, nor does it currently generate them for you. This must be performed manually by the user. See below for more information of the expected database files.

An example database sheet can look as follows, where 4 tools are being used, and `malt` and `kraken2` will be used against two databases each.

```console
tool,db_name,db_params,db_path
malt,malt85,-id 85,/<path>/<to>/malt/testdb-malt/
malt,malt95,-id 90,/<path>/<to>/malt/testdb-malt.tar.gz
kraken2,db1,,/<path>/<to>/kraken2/testdb-kraken2.tar.gz
kraken2,db2,--quick,/<path>/<to>/kraken2/testdb-kraken2.tar.gz
centrifuge,db1,,/<path>/<to>/centrifuge/minigut_cf.tar.gz
metaphlan3,db1,,/<path>/<to>/metaphlan3/metaphlan_database/
motus,db_mOTU,,/<path>/<to>/motus/motus_database/
```

Column specifications are as follows:

| Column      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tool`      | Taxonomic profiling tool (supported by nf-core/taxprofiler) that the database has been indexed for [required].                                                                                                                                                                                                                                                                                                                                                 |
| `db_name`   | A unique name of the particular database [required].                                                                                                                                                                                                                                                                                                                                                                                                           |
| `db_params` | Any parameters of the given taxonomic profiler that you wish to specify that the taxonomic profiling tool should use when profiling against this specific. Can be empty to use taxonomic profiler defaults. Must not be surrounded by quotes [required]. We generally do not recommend specifying parameters here that turn on/off saving of output files or specifying particular file extensions - this should be already addressed via pipeline parameters. |
| `db_path`   | Path to the database. Can either be a path to a directory containing the database index files or a `.tar.gz` file which contains the compressed database directory with the same name as the tar archive, minus `.tar.gz` [required].                                                                                                                                                                                                                          |

> 💡 You can also specify the same database directory/file twice (ensuring unique `db_name`s) and specify different parameters for each database to compare the effect of different parameters during profiling.

nf-core/taxprofiler will automatically decompress and extract any compressed archives for you.

Expected (uncompressed) database files for each tool are as follows:

- **MALT** output of `malt-build`. A directory containing:
  - `ref.idx`
  - `taxonomy.idx`
  - `taxonomy.map`
  - `index0.idx`
  - `table0.idx`
  - `table0.db`
  - `ref.inf`
  - `ref.db`
  - `taxonomy.tre`
- **Kraken2** output of `kraken2-build` command(s) A directory containing:
  - `opts.k2d`
  - `hash.k2d`
  - `taxo.k2d`
- **Centrifuge** output of `centrifuge-build`. A directory containing:
  - `<database_name>.<number>.cf`
  - `<database_name>.<number>.cf`
  - `<database_name>.<number>.cf`
  - `<database_name>.<number>.cf`
- **MetaPhlAn3** generated with `metaphlan --install` or downloaded from links on the [MetaPhlAn3 wiki](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0#customizing-the-database). A directory containing:
  - `mpa_v30_CHOCOPhlAn_201901.pkl`
  - `mpa_v30_CHOCOPhlAn_201901.pkl`
  - `mpa_v30_CHOCOPhlAn_201901.fasta`
  - `mpa_v30_CHOCOPhlAn_201901.3.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.4.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.1.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.2.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.rev.1.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.rev.2.bt2`
  - `mpa_latest`
- **Kaiju** output of `kaiju-makedb`. A directory containing:
  - `kaiju_db_*.fmi`
  - `nodes.dmp`
  - `names.dmp`
- **DIAMOND** output of `diamond makedb`. Note: requires building with taxonomy files
  to generate taxonomic profile. See [DIAMOND documentation](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#makedb-options). A file named:
  - `<database_name>.dmnd`
- **mOTUs** is composed of code and database together. The mOTUs tools
  [`downloadDB`](https://github.com/motu-tool/mOTUs/blob/master/motus/downloadDB.py)
  is used to prepare the mOTUs database and create a file with the version information.
  The database download step can be time consuming and the database will be consisting
  with same release version of the mOTUs tools. The database for same version tools
  can be thus reused for multiple runs. Users can download the database once using the script above and
  specify the path the database to the TSV table provided to `--databases`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/taxprofiler --input samplesheet.csv --databases databases.csv --outdir <OUTDIR> -profile docker --run_<TOOL1> --run_<TOOL2>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Preprocessing Steps

nf-core/taxprofiler offers four main preprocessing steps

- Read processing: adapter clipping and pair-merging.
- Complexity filtering: removal of low-sequence complexity reads.
- Host read-removal: removal of reads aligning to reference genome(s) of a host.
- Run merging: concatenation of multiple FASTQ chunks/sequencing runs/libraries of a sample.

#### Read Processing

Raw sequencing read processing in the form of adapter clipping and paired-end read merging can be activated via the `--perform_shortread_qc` or `--perform_longread_qc` flags.

It is highly recommended to run this on raw reads to remove artefacts from sequencing that can cause false positive identification of taxa (e.g. contaminated reference genomes) and/or skews in taxonomic abundance profiles.

There are currently two options for short-read preprocessing: `fastp` or `adapterremoval`.

For adapter clipping, you can either rely on tool default adapter sequences, or supply your own adapters (`--shortread_qc_adapter1` and `--shortread_qc_adapter2`)
By default, paired-end merging is not activated and paired-end profiling is performed where supported otherwise pairs will be independently profiled. If paired-end merging is activated you can also specify whether to exclude unmerged reads in the reads sent for profiling (`--shortread_qc_mergepairs` and `--shortread_qc_excludeunmerged`).
You can also turn off clipping and only perform paired-end merging, if requested. This can be useful when processing data downloaded from the ENA, SRA, or DDBJ (`--shortread_qc_skipadaptertrim`).
Both tools support length filtering of reads and can be tuned with `--shortread_qc_minlength`. Performing length filtering can be useful to remove short (often low sequencing complexity) sequences that result in unspecific classification and therefore slow down runtime during profiling, with minimal gain.

There is currently one option for long-read Oxford Nanopore processing: `porechop`.

For both short-read and long-read preprocessing, you can optionally save the resulting processed reads with `--save_preprocessed_reads`.

#### Complexity Filtering

Complexity filtering can be activated via the `--perform_shortread_complexityfilter` flag.

Complexity filtering is primarily a run-time optimisation step. It is not necessary for accurate taxonomic profiling, however it can speed up run-time of each tool by removing reads with low-diversity of nucleotides (e.g. with mono-nucleotide - `AAAAAAAA`, or di-nucleotide repeats `GAGAGAGAGAGAGAG`) that have a low-chance of giving an informative taxonomic ID as they can be associated with many different taxa. Removing these reads therefore saves computational time and resources.

There are currently three options for short-read complexity filtering: [`bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), [`prinseq++`](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus), and [`fastp`](https://github.com/OpenGene/fastp#low-complexity-filter).

The tools offer different algorithms and parameters for removing low complexity reads. We therefore recommend reviewing the pipeline's [parameter documentation](https://nf-co.re/taxprofiler/parameters) and the documentation of the tools (see links above) to decide on optimal methods and parameters for your dataset.

You can optionally save the FASTQ output of the run merging with the `--save_complexityfiltered_reads`. If running with `fastp`, complexity filtering happens inclusively within the earlier shortread preprocessing step. Therefore there will not be an independent pipeline step for complexity filtering, and no independent FASTQ file (i.e. `--save_complexityfiltered_reads` will be ignored) - your complexity filtered reads will also be in the `fastp/` folder in the same file(s) as the preprocessed read.

#### Host Removal

Removal of possible-host reads from FASTQ files prior profiling can be activated with `--perform_shortread_hostremoval` or `--perform_longread_hostremoval`.

Similarly to complexity filtering, host-removal can be useful for runtime optimisation and reduction in misclassified reads. It is not always necessary to report classification of reads from a host when you already know the host of the sample, therefore you can gain a run-time and computational advantage by removing these prior typically resource-heavy profiling with more efficient methods. Furthermore, particularly with human samples, you can reduce the number of false positives during profiling that occur due to host-sequence contamination in reference genomes on public databases.

nf-core/taxprofiler currently offers host-removal via alignment against a reference genome with Bowtie2, and the use of the unaligned reads for downstream profiling.

You can supply your reference genome in FASTA format with `--hostremoval_reference`. You can also optionally supply a directory containing pre-indexed Bowtie2 index files with `--shortread_hostremoval_index` or a minimap2 `.mmi` file for `--longread_hostremoval_index`, however nf-core/taxprofiler will generate these for you if necessary. Pre-supplying the index directory or files can greatly speed up the process, and these can be re-used.

> 💡 If you have multiple taxa or sequences you wish to remove (e.g., the host genome and then also PhiX - common quality-control reagent during sequencing) you can simply concatenate the FASTAs of each taxa or sequences into a single reference file.

#### Run Merging

For samples that may have been sequenced over multiple runs, or for FASTQ files split into multiple chunks, you can activate the ability to merge across all runs or chunks with `--perform_runmerging`.

For more information how to set up your input samplesheet, see [Multiple runs of the same sample](#multiple-runs-of-the-same-sample).

Activating this functionality will concatenate the FASTQ files with the same sample name _after_ the optional preprocessing steps and _before_ profiling. Note that libraries with runs of different pairing types will **not** be merged and this will be indicated on output files with a `_se` or `_pe` suffix to the sample name accordingly.

You can optionally save the FASTQ output of the run merging with the `--save_runmerged_reads`.

##### Profiling

###### MALT

nf-core/taxprofiler uses MALT 0.4.1, which is a compatively old version. However it has been found that the most recent version of MALT (0.5.\*), at the time of writing, is broken. [The the LCA step appears not to be executed](http://megan.informatik.uni-tuebingen.de/t/lca-placement-failure-with-malt-v-0-5-2-and-0-5-3/1996/3), pushing all hits to the leaves of the taxonomy. However, if you need to use a more recent taxonomy map file with your databases, the output of `malt-build` from MALT 0.5.3 should be still be compatible with `malt-run` of 0.4.1.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/taxprofiler
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/taxprofiler releases page](https://github.com/nf-core/taxprofiler/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

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
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

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

## Troubleshooting and FAQs

### I get a warning during centrifuge_kreport process with exit status 255.

When a sample has insufficient hits for abundance estimation, the resulting `report.txt` file will be empty.

When trying to convert this to a kraken-style report, the conversion tool will exit with a status code `255`, and provide a `WARN`.

This is **not** an error nor a failure of the pipeline, just your sample has no hits to the provided database when using centrifuge.
