# nf-core/taxprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [date]

### `Added`

- [#298](https://github.com/nf-core/taxprofiler/pull/298) **New classifier** [ganon](https://pirovc.github.io/ganon/) (added by @jfy133)
- [#276](https://github.com/nf-core/taxprofiler/pull/276) Implemented batching in the KrakenUniq samples processing (added by @Midnighter)
- [#272](https://github.com/nf-core/taxprofiler/pull/272) Add saving of final 'analysis-ready-reads' to dedicated directory (❤️ to @alexhbnr for request, added by @jfy133)
- [#303](https://github.com/nf-core/taxprofiler/pull/303) Add support for taxpasta profile standardisation in single sample pipeline runs (❤️ to @artur-matysik for request, added by @jfy133)
- [#308](https://github.com/nf-core/taxprofiler/pull/308) Add citations and bibliographic information to the MultiQC methods text of tools used in a given pipeline run (added by @jfy133)
- [#315](https://github.com/nf-core/taxprofiler/pull/315) Updated to nf-core pipeline template v2.9 (added by @sofstam & @jfy133)
- [#319](https://github.com/nf-core/taxprofiler/pull/319) Added support for virus hit expansion in Kaiju (❤️ to @dnlrxn for requesting, added by @jfy133)
- [#323](https://github.com/nf-core/taxprofiler/pull/323) Add ability to skip sequencing quality control tools (❤️ to @vinisalazar for requesting, added by @jfy133)
- [#318](https://github.com/nf-core/taxprofiler/pull/318) Added the profiler MetaPhlAn4 and removed MetaPhlAn3 (added by @LilyAnderssonLee)

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
- [#307](https://github.com/nf-core/taxprofiler/pull/307) Fix databases being sometimes associated with the wrong tool (e.g. Kaiju) (fix by @jfy133)
- [#313](https://github.com/nf-core/taxprofiler/pull/304) Fix pipeline not providing error when database sheet does not have a header (♥ to @noah472 for reporting, fix by @jfy133)
- [#330](https://github.com/nf-core/taxprofiler/pull/330) Added better tagging to allow disambiguation of Kraken2 steps of Kraken2 vs Bracken (♥ to @MajoroMask for requesting, added by @jfy133)
- [#334](https://github.com/nf-core/taxprofiler/pull/334) Increase the memory of the FALCO process to 4GB (fix by @LilyAnderssonLee)
- [#332](https://github.com/nf-core/taxprofiler/pull/332) Improved meta map stability for more robust pipeline resuming (fix by @jfy133)
- [#338](https://github.com/nf-core/taxprofiler/pull/338) Fixed wrong file 'out' file going to `centrifuge kreport` module (♥ to @LilyAnderssonLee for reporting, fix by @jfy133)
- [#342](https://github.com/nf-core/taxprofiler/pull/342) Fixed docs/usage to correctly list the required database files for Bracken and tips to obtain Kraken2 databases (fix by @husensofteng)

### `Dependencies`

| Tool      | Previous version | New version |
| --------- | ---------------- | ----------- |
| MultiQC   | 1.13             | 1.14        |
| taxpasta  | 0.2.3            | 0.4.1       |
| MetaPhlAn | 3.0.12           | 4.0.6       |

### `Deprecated`

- [#338](https://github.com/nf-core/taxprofiler/pull/338) Updated Centrifuge module to not generate (undocumented) SAM alignments by default if --save_centrifuge_reads supplied, to due to a Centrifuge bug modifying profile header. SAM alignments can still be generated if `--out-fmt` supplied in `database.csv` (♥ to @LilyAnderssonLee for reporting, fix by @jfy133)

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
