# Troubleshooting and FAQs

## Which sequencing read QC tool should I use?

nf-core/taxprofiler offers two options for quality control reporting of sequencing reads, [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Falco](https://github.com/smithlabcode/falco).

Falco is a recent, more performant C++, implementation of the Java-based FastQC, while producing the same output files.
Falco is starting to add additional features, but is not as widely used as FastQC.

We offer FastQC as the gold-standard and most used option within the genomics community.
Falco is included as a faster alternative for users who wish to try it out or need to save time and resources.

## Which short read-adapter clipping/read merging tool should I use?
`
nf-core/taxprofiler offers two options for short-read adapter clipping and read merging, [fastp](https://github.com/OpenGene/fastp) and [AdapterRemoval](https://github.com/MikkelSchubert/AdapterRemoval) (v2).

The two tools are both well-built tools, and offered based on popularity by different communities.
For example, AdapterRemoval (v2) is often used in ancient DNA studies due to additional functionality for additional barcode data, while fastp is commonly used in modern metagenomics studies with additional general sequencing data cleanup functionality.

Each have their own advantages and disadvantages, and the choice of which to use is up to the user.
You should evaluate both tools, through literature review (such as benchmarking studies), and/or by testing on your own data.

## Which long-read trimming tool should I use?

nf-core/taxprofiler offers two options for Oxford Nanopore Technologies (ONT) data long-read adapter clipping and read trimming, [Porechop](https://github.com/rrwick/Porechop) and [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI).

The two tools are very similar, with the main difference being that Porechop_ABI is a more recent fork of the original Porechop, with the same functionality but without the need to pre-define adapter sequences.

Selection of the tool should be based on reading the documentation and comparing the features to the structure of your own data.

## Which long-read length filtering tool should I use?

nf-core/taxprofiler offers two options for long-read length filtering, [nanoq](https://github.com/esteinig/nanoq) and [Filtlong](https://github.com/rrwick/Filtlong).

The two tools are designed to perform the same task of filtering ONT data long reads.

nanoq was developed as a more actively developed and a more performant alternative to Filtlong. It is written in Rust, with performance reported to be much faster and more efficient than the older Filtlong tool.
However, Filtlong has started to be slightly more active again.

Each have their own advantages and disadvantages, and the choice of which to use is up to the user.
You should evaluate both tools, through literature review (such as benchmarking studies), and/or by testing on your own data.

## Which short-read complexity filtering tool should I use?

nf-core/taxprofiler offers three options for short-read complexity filtering, [PRINSEQ++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus), [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/), and [fastp](https://github.com/OpenGene/fastp).

PRINSEQ++ is a recent and more efficient C++ version of the original PRINSEQ tool that was written in Perl.
BBDuk is a well-established and widely used Java based toolkit with a large number of features and a strong user base.
fastp is a general purpose sequencing QC tool that provides complexity filtering as an integrated step of read-merging and adapter trimming.

The three tools all perform well, and are offered as each provide different algorithms for complexity filtering.

- **PRINSEQ++**: filter based on an estimate of sequence randomness (entropy) using a [DUST score](https://pubmed.ncbi.nlm.nih.gov/16796549/).
- **BBDuk**: filter based on k-mer frequency
- **fastp**: also includes a complexity filtering step integrated into its own adapter-trimming steps

Each have their own advantages and disadvantages, with different tuning parameters and the choice of which to use is up to the user.
You should evaluate all three tools, through literature review (such as benchmarking studies), and/or by testing on your own data.

## Which short-read host-removal tool should I use?

nf-core/taxprofiler offers three different tools for host-removal of short-reads.

- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [Hostile](https://github.com/bede/hostile)
- [Deacon](https://github.com/bede/deacon)

Hostile and Bowtie2 provide alignment-based host-removal, while Deacon provides a k-mer based approach.
If you wish to use an alignment approach for the removal of host DNA, and are working with Human-associated or common-model organism samples, we recommend using Hostile over Bowtie2.
Hostile used Bowtie2 internally but as it provides optimised alignment presets as well as superior reference genomes (such as the T2T CHM13v2.0 human genome).

Deacon is a newer k-mer based approach that is reportedly much faster and efficient than alignment-based methods, with the efficiency improvements allowing the use of more inclusive pan-genomes for host-removal.

We recommend evaluating between Hostile and Deacon through literature review (such as benchmarking studies), and/or by testing on your own data.

## I get a warning during CENTRIFUGE_KREPORT or KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE process with exit status 255

When a sample has insufficient hits for abundance estimation, the resulting `report.txt` file will be empty.

When trying to convert this to a kraken-style report or merging together, the tools will exit with a status code `255`, and provide a `WARN`.

This is _not_ an error nor a failure of the pipeline, just your sample has no hits to the provided database when using Centrifuge.

## Why does any error or failed process from BRACKEN_BRACKEN get ignored?

If Kraken2 doesn't classify any reads (100% unclassified) Bracken will fail on that input.
Having no reads taxonomically classified can be a reasonable outcome for some samples, and we don't want to treat this as an error.

Therefore to allow the remainder of a given run to complete if a single Bracken task fails, [we set the Nextflow `errorStrategy` to `ignore`](https://github.com/nf-core/taxprofiler/blob/5d3ee5513a84f92773c8376c55b5f4da39835307/conf/base.config#L61-L63).
You will still get a warning in the Nextflow console output and log that a Bracken task failed.

If you want to change this behaviour, you can override the `errorStrategy` in your own Nextflow configuration file.

For example:

```nextflow
process {
    withName: BRACKEN_BRACKEN {
        errorStrategy = 'retry'
    }
}
```

Other `errorStrategy` options can be found on the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#errorstrategy).
