# nf-core/taxprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [date]

### `Added`

- [#276](https://github.com/nf-core/taxprofiler/pull/276) Implemented batching in the KrakenUniq samples processing. (added by @Midnighter)
- [#272](https://github.com/nf-core/taxprofiler/pull/272) - Add saving of final 'analysis-ready-reads' to dedicated directory. (❤️ to @alexhbnr for reporting, added by @jfy133)

### `Fixed`

- [#271](https://github.com/nf-core/taxprofiler/pull/271/files) Improved standardised table generation documentation nd mOTUs manual database download tutorial (♥ to @prototaxites for reporting, fix by @jfy133)
- [#269](https://github.com/nf-core/taxprofiler/pull/269/files) Reduced output files in AWS full test output due to very large files
- [#270](https://github.com/nf-core/taxprofiler/pull/270/files) Fixed warning for host removal index parameter, and improved index checks (♥ to @prototaxites for reporting, fix by @jfy133)
- [#274](https://github.com/nf-core/taxprofiler/pull/274/files) Substituted the samtools/bam2fq module with samtools/fastq module (fix by @sofstam)
- [#275](https://github.com/nf-core/taxprofiler/pull/275/files) Replaced function used for error reporting to more Nextflow friendly method (fix by @jfy133)
- [#285](https://github.com/nf-core/taxprofiler/pull/285/files) Fixed overly large log files in Kraken2 output (♥ to @prototaxites for reporting, fix by @Midnighter & @jfy133)
- [#286](https://github.com/nf-core/taxprofiler/pull/286/files) Runtime optimisation of MultiQC step via improved log file processing (fix by @Midnighter & @jfy133)
- [#289](https://github.com/nf-core/taxprofiler/pull/289/files) Pipeline updated to nf-core template 2.8 (fix by @Midnighter & @jfy133)
- [#290](https://github.com/nf-core/taxprofiler/pull/286/files) Minor database input documentation improvements (♥ to @alneberg for reporting, fix by @jfy133)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| MultiQC | 1.13             | 1.14        |

### `Deprecated`

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
