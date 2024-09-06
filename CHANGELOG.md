# nf-core/taxprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2dev - Bouncy Basenji [unreleased]

### `Added`

- [#417](https://github.com/nf-core/taxprofiler/pull/417) - Added reference-free metagenome estimation with Nonpareil (added by @jfy133)
- [#466](https://github.com/nf-core/taxprofiler/pull/466) - Input database sheets now require a `db_type` column to distinguish between short- and long-read databases (added by @LilyAnderssonLee)
- [#505](https://github.com/nf-core/taxprofiler/pull/505) - Add small files to the file `tower.yml` (added by @LilyAnderssonLee)
- [#508](https://github.com/nf-core/taxprofiler/pull/508) - Add `nanoq` as a filtering tool for nanopore reads (added by @LilyAnderssonLee)
- [#511](https://github.com/nf-core/taxprofiler/pull/511) - Add `porechop_abi` as an alternative adapter removal tool for long reads nanopore data (added by @LilyAnderssonLee)
- [#512](https://github.com/nf-core/taxprofiler/pull/512) - Update all tools to the latest version and include nf-test (Updated by @LilyAnderssonLee & @jfy133)

### `Fixed`

- [#518](https://github.com/nf-core/taxprofiler/pull/518) Fixed a bug where Oxford Nanopore FASTA input files would not be processed (❤️ to @ikarls for reporting, fixed by @jfy133)

### `Dependencies`

| Tool                     | Previous version | New version |
| ------------------------ | ---------------- | ----------- |
| bbmap                    | 39.01            | 39.06       |
| bowtie2                  | 2.4.4            | 2.5.2       |
| bracken                  | 2.7              | 2.9         |
| cat/fastq                | 8.30             |
| diamond                  | 2.0.15           | 2.1.8       |
| ganon                    | 1.5.1            | 2.0.0       |
| kraken2                  | 2.1.2            | 2.1.3       |
| krona                    | 2.8              | 2.8.1       |
| megan                    | 6.24.20          | 6.25.9      |
| metaphlan                | 4.0.6            | 4.1.1       |
| minimap2                 | 2.24             | 2.28        |
| motus/profile            | 3.0.3            | 3.1.0       |
| multiqc                  | 1.21             | 1.24.1      |
| nanoq                    |                  | 0.10.0      |
| samtools                 | 1.17             | 1.20        |
| untar                    | 4.7              | 4.8         |

### `Deprecated`

## v1.1.8 - Augmented Akita Patch [2024-06-20]

### `Added`

- [#487](https://github.com/nf-core/taxprofiler/pull/487) Updated to nf-core pipeline template v2.14.1 (added by @jfy133)

### `Fixed`

- [#484](https://github.com/nf-core/taxprofiler/pull/484) Improved input validation to immediately fail if run accession IDs within a given sample ID are not unique (❤️ to @sofstam for reporting, fixed by @jfy133)
- [#491](https://github.com/nf-core/taxprofiler/pull/491) Added flag to publish intermediate bracken files (❤️ to @ewissel for reporting, fixed by @sofstam and @jfy133)
- [#489](https://github.com/nf-core/taxprofiler/pull/489) Fix KrakenUniq classified reads output format mismatch (❤️ to @SannaAb for reporting, fixed by @jfy133)
- [#495](https://github.com/nf-core/taxprofiler/pull/495) Stop TAXPASTA failures when profiles do not have exact compositionality (fixes by @Midnighter, @jfy133)

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| KMCP     | 0.9.1            | 0.9.4       |
| TAXPASTA | 0.6.1            | 0.7.0       |

### `Deprecated`

- [#492](https://github.com/nf-core/taxprofiler/pull/492) Removed `--kmcp_mode` parameter from KMCP to allow per database specification by setting in db_params in database sheet (fixed by @jfy133)

## v1.1.7 - Augmented Akita Patch [2024-04-25]

### `Added`

- [#477](https://github.com/nf-core/taxprofiler/pull/477) Provide more emphasis and links to tutorials on how to retrieve and supply reference databases (❤️ to @vmkalbskopf for reporting, added by @jfy133)

### `Fixed`

- [#476](https://github.com/nf-core/taxprofiler/pull/476/) Fixed bug in validating Bracken/Kraken/KMCP split database parameters (fixed by @LilyAnderssonLee)

### `Dependencies`

### `Deprecated`

## v1.1.6 - Augmented Akita Patch [2024-04-16]

### `Added`

- [#454](https://github.com/nf-core/taxprofiler/pull/454) Updated to nf-core pipeline template v2.13.1 (added by @LilyAnderssonLee & @sofstam)
- [#461](https://github.com/nf-core/taxprofiler/pull/461) Turned on 'strict' Nextflow evaluation runs (added by @jfy133)
- [#461](https://github.com/nf-core/taxprofiler/pull/461) Optimised database compression so each compressed input database is untarred once, and shared amongst each run with different parameters (added by @jfy133)
- [#461](https://github.com/nf-core/taxprofiler/pull/461) Added new parameter to optionally save uncompressed databases (added by @jfy133)
- [#471](https://github.com/nf-core/taxprofiler/pull/471) Remove `-stub` run in the `download_pipeline.yml` because the pipeline does not support stub runs on dev (fixed by @LilyAnderssonLee)

### `Fixed`

- [#336](https://github.com/nf-core/taxprofiler/issues/336) Replace samplesheet check with nf-validation for both sample and database input sheets (fix by @LilyAnderssonLee)
- [#460](https://github.com/nf-core/taxprofiler/issues/460) corrected the channel transformations to combine Kaiju and mOTUs reports with their reference databases (fix by @Midnighter)

### `Dependencies`

### `Deprecated`

## v1.1.5 - Augmented Akita Patch [2024-02-08]

### `Added`

- [#439](https://github.com/nf-core/taxprofiler/pull/439) Read deduplication with fastp (added by @maxibor)
- [#440](https://github.com/nf-core/taxprofiler/pull/440) Include mention of pre-built kaiju databases in tutorial.md (added by @Joon-Klaps)
- [#442](https://github.com/nf-core/taxprofiler/pull/442) Updated to nf-core pipeline template v2.12 (added by @sofstam)

### `Fixed`

- [#444](https://github.com/nf-core/taxprofiler/pull/444) Centrifuge now uses dedicated tmp directory to hopefully prevent mkfifo clashes (❤️ to @erinyoung for reporting, fix by @jfy133)

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| Centrifuge | 1.0.4_beta       | 1.0.4.1     |

### `Deprecated`

## v1.1.4 - Augmented Akita Patch [2024-01-24]

### `Added`

### `Fixed`

- [#431](https://github.com/nf-core/modules/pull/4781#event-11555493525) Updated kaiju2table module to report taxon names (fix by @Joon-Klaps)
- [#430](https://github.com/nf-core/taxprofiler/pull/430) Fix the fastq output in the module LONGREAD_HOSTREMOVAL. (fix by @LilyAnderssonLee)

### `Dependencies`

| Tool  | Previous version | New version |
| ----- | ---------------- | ----------- |
| kaiju | 1.8.2            | 1.10.0      |

### `Deprecated`

## v1.1.3 - Augmented Akita Patch [2024-01-12]

### `Added`

- [#424](https://github.com/nf-core/taxprofiler/pull/424) Updated to nf-core pipeline template v2.11.1 (added by @LilyAnderssonLee & @sofstam)

### `Fixed`

- [#419](https://github.com/nf-core/taxprofiler/pull/419) Added improved syntax highlighting for tables in documentation (fix by @mashehu)
- [#421](https://github.com/nf-core/taxprofiler/pull/421) Updated the krakenuniq/preloadedkrakenuniq module that contained a fix for saving the output reads (❤️ to @SannaAb for reporting, fix by @Midnighter)
- [#427](https://github.com/nf-core/taxprofiler/pull/427) Fixed preprint information in the recommended methods text (fix by @jfy133)

### `Dependencies`

| Tool          | Previous version | New version |
| ------------- | ---------------- | ----------- |
| multiqc       | 1.15             | 1.19        |
| fastqc        | 11.9             | 12.1        |
| nf-validation | unpinned         | 1.1.3       |

## v1.1.2 - Augmented Akita Patch [2023-10-27]

### `Added`

- [#408](https://github.com/nf-core/taxprofiler/pull/408) Added preprint citation information to README and manifest (added by @jfy133)

### `Fixed`

- [#405](https://github.com/nf-core/taxprofiler/pull/405) Fix database to tool mismatching in KAIJU2KRONA input (❤️ to @MajoroMask for reporting, fix by @jfy133)
- [#406](https://github.com/nf-core/taxprofiler/pull/406) Fix overwriting of bracken-derived kraken2 outputs when the database name is shared between Bracken/Kraken2. (❤️ to @MajoroMask for reporting, fix by @jfy133)
- [#409](https://github.com/nf-core/taxprofiler/pull/409) Fix a NullPointerException error occurring occasionally in older version of MEGAN's rma2info (❤️ to @MajoroMask for reporting, fix by @jfy133)

### `Dependencies`

| Tool           | Previous version | New version |
| -------------- | ---------------- | ----------- |
| megan/rma2info | 6.21.7           | 6.24.20     |

### `Deprecated`

## v1.1.1 - Augmented Akita Patch [2023-10-11]

### `Added`

- [#379](https://github.com/nf-core/taxprofiler/pull/379) Added support for previously missing Bracken-corrected Kraken2 report as output (added by @hkaspersen & @jfy133 )
- [#380](https://github.com/nf-core/taxprofiler/pull/380) Updated to nf-core pipeline template v2.10 (added by @LilyAnderssonLee & @sofstam)
- [#393](https://github.com/nf-core/taxprofiler/pull/383) Add validation check for a taxpasta taxonomy directory if --taxpasta*add*\* parameters requested (♥️ to @alimalrashed for reporting, added by @jfy133)

### `Fixed`

- [#383](https://github.com/nf-core/taxprofiler/pull/383) Update the module of KrakenUniq to the latest to account for edge case bugs where FASTQ input was mis-detected as wrong format (❤️ to @asafpr for reporting and solution, fixed by @LilyAnderssonLee)
- [#392](https://github.com/nf-core/taxprofiler/pull/392) Update the module of Taxpasta to support adding taxa information to results (❤️ to @SannaAb for reporting, fixed by @Midnighter)

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| KrakenUniq | 1.0.2            | 1.0.4       |
| taxpasta   | 0.6.0            | 0.6.1       |

### `Deprecated`

## v1.1.0 - Augmented Akita [2023-09-19]

### `Added`

- [#298](https://github.com/nf-core/taxprofiler/pull/298) **New classifier** [ganon](https://pirovc.github.io/ganon/) (added by @jfy133)
- [#312](https://github.com/nf-core/taxprofiler/pull/312) **New classifier** [KMCP](https://github.com/shenwei356/kmcp) (added by @sofstam)
- [#318](https://github.com/nf-core/taxprofiler/pull/318) **New classifier** [MetaPhlAn4](https://github.com/biobakery/MetaPhlAn) (MetaPhlAn3 support remains) (added by @LilyAnderssonLee)
- [#276](https://github.com/nf-core/taxprofiler/pull/276) Implemented batching in the KrakenUniq samples processing (added by @Midnighter)
- [#272](https://github.com/nf-core/taxprofiler/pull/272) Add saving of final 'analysis-ready-reads' to dedicated directory (❤️ to @alexhbnr for request, added by @jfy133)
- [#303](https://github.com/nf-core/taxprofiler/pull/303) Add support for taxpasta profile standardisation in single sample pipeline runs (❤️ to @artur-matysik for request, added by @jfy133)
- [#308](https://github.com/nf-core/taxprofiler/pull/308) Add citations and bibliographic information to the MultiQC methods text of tools used in a given pipeline run (added by @jfy133)
- [#315](https://github.com/nf-core/taxprofiler/pull/315) Updated to nf-core pipeline template v2.9 (added by @sofstam & @jfy133)
- [#321](https://github.com/nf-core/taxprofiler/pull/321) Added support for virus hit expansion in Kaiju (❤️ to @dnlrxn for requesting, added by @jfy133)
- [#325](https://github.com/nf-core/taxprofiler/pull/325) Add ability to skip sequencing quality control tools (❤️ to @vinisalazar for requesting, added by @jfy133)
- [#345](https://github.com/nf-core/taxprofiler/pull/345) Add simple tutorial to explain how to get up and running with an nf-core/taxprofiler run (added by @jfy133)
- [#355](https://github.com/nf-core/taxprofiler/pull/355) Add support for TAXPASTA's `--add-rank-lineage` to output (❤️ to @MajoroMask for request, added by @Midnighter, @sofstam, @jfy133)
- [#368](https://github.com/nf-core/taxprofiler/pull/368/) Add the ability to ignore profile errors caused by empty profiles and other validation errors when merging multiple profiles using TAXPASTA (added by @Midnighter and @LilyAnderssonLee)

### `Fixed`

- [#271](https://github.com/nf-core/taxprofiler/pull/271) Improved standardised table generation documentation for mOTUs manual database download tutorial (♥ to @prototaxites for reporting, fix by @jfy133)
- [#269](https://github.com/nf-core/taxprofiler/pull/269) Reduced output files in AWS full test output due to very large files (fix by @jfy133)
- [#270](https://github.com/nf-core/taxprofiler/pull/270) Fixed warning for host removal index parameter, and improved index checks (♥ to @prototaxites for reporting, fix by @jfy133)
- [#274](https://github.com/nf-core/taxprofiler/pull/274) Substituted the samtools/bam2fq module with samtools/fastq module (fix by @sofstam)
- [#275](https://github.com/nf-core/taxprofiler/pull/275) Replaced function used for error reporting to more Nextflow friendly method (fix by @jfy133)
- [#285](https://github.com/nf-core/taxprofiler/pull/285) Fixed overly large log files in Kraken2 output (♥ to @prototaxites for reporting, fix by @Midnighter & @jfy133)
- [#286](https://github.com/nf-core/taxprofiler/pull/286) Runtime optimisation of MultiQC step via improved log file processing (fix by @Midnighter & @jfy133)
- [#289](https://github.com/nf-core/taxprofiler/pull/289) Pipeline updated to nf-core template 2.8 (fix by @Midnighter & @jfy133)
- [#290](https://github.com/nf-core/taxprofiler/pull/290) Minor database input documentation improvements (♥ to @alneberg for reporting, fix by @jfy133)
- [#305](https://github.com/nf-core/taxprofiler/pull/305) Fix docker/podman registry definition for tower compatibility (fix by @adamrtalbot, @jfy133)
- [#304](https://github.com/nf-core/taxprofiler/pull/304) Correct mistake in kaiju2table documentation, only single rank can be supplied (♥ to @artur-matysik for reporting, fix by @jfy133)
- [#307](https://github.com/nf-core/taxprofiler/pull/307) Fix databases being sometimes associated with the wrong tool (e.g. Kaiju) (fix by @jfy133, @Midnighter and @LilyAnderssonLee)
- [#313](https://github.com/nf-core/taxprofiler/pull/313) Fix pipeline not providing error when database sheet does not have a header (♥ to @noah472 for reporting, fix by @jfy133)
- [#330](https://github.com/nf-core/taxprofiler/pull/330) Added better tagging to allow disambiguation of Kraken2 steps of Kraken2 vs Bracken (♥ to @MajoroMask for requesting, added by @jfy133)
- [#334](https://github.com/nf-core/taxprofiler/pull/334) Increase the memory of the FALCO process to 4GB (fix by @LilyAnderssonLee)
- [#332](https://github.com/nf-core/taxprofiler/pull/332) Improved meta map stability for more robust pipeline resuming (fix by @jfy133)
- [#338](https://github.com/nf-core/taxprofiler/pull/338) Fixed wrong file 'out' file going to `centrifuge kreport` module (♥ to @LilyAnderssonLee for reporting, fix by @jfy133)
- [#342](https://github.com/nf-core/taxprofiler/pull/342) Fixed docs/usage to correctly list the required database files for Bracken and tips to obtain Kraken2 databases (fix by @husensofteng)
- [#350](https://github.com/nf-core/taxprofiler/pull/350) Reorganize the CI tests into separate profiles in preparation for implementation of nf-test (fix by @LilyAnderssonLee)
- [#364](https://github.com/nf-core/taxprofiler/pull/364) Add autoMounts to apptainer profile in nextflow.config (♥ to @hkaspersen for reporting, fix by @LilyAnderssonLee)
- [#372](https://github.com/nf-core/taxprofiler/pull/372) Update modules to use quay.io nf-core mirrored containers (♥ to @maxulysse for pointing out, fix by @LilyAnderssonLee and @jfy133)

### `Dependencies`

| Tool      | Previous version | New version |
| --------- | ---------------- | ----------- |
| MultiQC   | 1.13             | 1.15        |
| TAXPASTA  | 0.2.3            | 0.6.0       |
| MetaPhlAn | 3.0.12           | 4.0.6       |
| fastp     | 0.23.2           | 0.23.4      |
| samtools  | 1.16.1           | 1.17        |

### `Deprecated`

- [#338](https://github.com/nf-core/taxprofiler/pull/338) Updated Centrifuge module to not generate (undocumented) SAM alignments by default if --save_centrifuge_reads supplied, due to a Centrifuge bug modifying profile header. SAM alignments can still be generated if `--out-fmt` supplied in `database.csv` (♥ to @LilyAnderssonLee for reporting, fix by @jfy133)

## v1.0.1 - Dodgy Dachshund Patch [2023-05-15]

### `Added`

### `Fixed`

- [#291](https://github.com/nf-core/taxprofiler/pull/291) - Fix Taxpasta not receiving taxonomy directory (❤️ to @SannaAb for reporting, fix by @jfy133)

## v1.0.0 - Dodgy Dachshund [2023-03-13]

Initial release of nf-core/taxprofiler, created with the [nf-core](https://nf-co.re/) template.

- Add read quality control (sequencing QC, adapter removal and merging)
- Add read complexity filtering
- Add host-reads removal step
- Add run merging
- Add taxonomic classification
- Add taxon table standardisation
- Add post-classification visualisation

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
