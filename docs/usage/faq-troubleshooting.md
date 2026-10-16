# Troubleshooting and FAQs

## Which classifier(s) should I use?

Choosing a classifier is not easy!
There is a reason why nf-core/taxprofiler is designed to offer multiple classifiers.
There are main pros- and cons- to each, and the best choice will depend on your research question, the information you need, and the trade-offs needed regarding resources (amongst others).

Here we provide a brief summary of the classifiers available in nf-core/taxprofiler.
With this information we hope to help you make an informed decision on which classifier(s) to use for your particular use case analysis.

| Classifier | Sequence read type | Database type         | Primary algorithm | Unique characteristic                                                                          | Common applications                                           |
| ---------- | ------------------ | --------------------- | ----------------- | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| Bracken    | short              | Nucleotide            | k-mer             | Redistribution of ambiguously classified reads by Kraken2 using a Bayesian model               | metagenomics/abundance estimation for closely related species |
| Centrifuge | short              | Nucleotide            | FM-index          | Does not require a lot of RAM like the other k-mer based classifiers                           | metagenomics, pathogen detection                              |
| DIAMOND    | short/long         | Amino acid            | alignment based   | High sensitivity and it runs much faster than the traditional BLAST                            | metagenomics, pathogen detection                              |
| ganon      | short              | Nucleotide            | k-mer             | Can update the databases with new sequencing without needing to rebuild the index from scratch | metagenomics                                                  |
| Kaiju      | short              | Amino acid            | alignment based   | Fast taxonomic classification at the protein level                                             | metagenomics, microbial ecology                               |
| KMCP       | short              | Nucleotide            | k-mer             | Reduces the high false-positive rates compared to other k-mer based classifiers                | pathogen detection                                            |
| Kraken2    | short              | Nucleotide            | k-mer             | Ultrafast speed at a minimal memory usage                                                      | pathogen detection                                            |
| KrakenUniq | short              | Nucleotide            | k-mer             | Reduced false-positive taxonomic classifications                                               | pathogen detection                                            |
| MALT       | short              | Nucleotide/Amino acid | alignment based   | It can provide functional classification along taxonomic classification                        | ancient metagenomics, microbiome profiling                    |
| melon      | long               | Nucleotide            | alignment based   | Two steps classification to handle long-read sequencing data                                   | microbiome profiling                                          |
| MetaPhlAn  | short              | Nucleotide            | alignment based   | Uses clade-specific marker gene                                                                | microbiome profiling                                          |
| MetaCache  | short/long         | Nucleotide            | k-mer             | Uses subset of k-mers using minhashing method and reduces memory requirements                  | metagenomics, pathogen detection, microbial ecology           |
| mOTUs      | short/long         | Nucleotide            | alignment based   | Uses universal marker genes and it does not rely on complete reference genomes                 | microbiome profiling                                          |
| sylph      | short/long         | Nucleotide            | k-mer             | Calculates Average Nucleotide Identity (ANI) and includes rapid processing time                | metagenomics, pathogen detection                              |

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
