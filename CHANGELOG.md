# nf-core/taxprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [date]

### `Added`

- [#276](https://github.com/nf-core/taxprofiler/pull/276) Implemented batching in the KrakenUniq samples processing. (added by @Midnighter)
- [#272](https://github.com/nf-core/taxprofiler/pull/272) Add saving of final 'analysis-ready-reads' to dedicated directory. (❤️ to @alexhbnr for reporting, added by @jfy133)

### `Fixed`

- [#271](https://github.com/nf-core/taxprofiler/pull/271) Improved standardised table generation documentation for mOTUs manual database download tutorial (♥ to @prototaxites for reporting, fix by @jfy133)
- [#269](https://github.com/nf-core/taxprofiler/pull/269) Reduced output files in AWS full test output due to very large files
- [#270](https://github.com/nf-core/taxprofiler/pull/270) Fixed warning for host removal index parameter, and improved index checks (♥ to @prototaxites for reporting, fix by @jfy133)
- [#274](https://github.com/nf-core/taxprofiler/pull/274) Substituted the samtools/bam2fq module with samtools/fastq module (fix by @sofstam)
- [#275](https://github.com/nf-core/taxprofiler/pull/275) Replaced function used for error reporting to more Nextflow friendly method (fix by @jfy133)
- [#285](https://github.com/nf-core/taxprofiler/pull/285) Fixed overly large log files in Kraken2 output (♥ to @prototaxites for reporting, fix by @Midnighter & @jfy133)
- [#286](https://github.com/nf-core/taxprofiler/pull/286) Runtime optimisation of MultiQC step via improved log file processing (fix by @Midnighter & @jfy133)
- [#289](https://github.com/nf-core/taxprofiler/pull/289) Pipeline updated to nf-core template 2.8 (fix by @Midnighter & @jfy133)
- [#290](https://github.com/nf-core/taxprofiler/pull/290) Minor database input documentation improvements (♥ to @alneberg for reporting, fix by @jfy133)
- [#305](https://github.com/nf-core/taxprofiler/pull/305) Fix docker/podman registry definition for tower compatibility (fix by @adamrtalbot, @jfy133)
- [#304](https://github.com/nf-core/taxprofiler/pull/304) Correct mistake in kaiju2table documentation, only single rank can be supplied (♥ to @artur-matysik for reporting, fix by @jfy133)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| MultiQC | 1.13             | 1.14        |

### `Deprecated`

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
