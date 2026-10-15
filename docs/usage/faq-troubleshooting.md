# Troubleshooting and FAQs

## Which classifier(s) should I use?

Choosing a classifier is not easy!
There is a reason why nf-core/taxprofiler is designed to offer multiple classifiers.
There are main pros- and cons- to each, and the best choice will depend on your research question, the information you need, and the trade-offs needed regarding resources (amongst others).

Here we provide a brief summary of the classifiers available in nf-core/taxprofiler.
With this information we hope to help you make an informed decision on which classifier(s) to use for your particular use case analysis.

| Classifier | Sequence read type | Database type          | Primary algorithm | Unique characteristic |
| ---------- | ------------------ | ---------------------- | ----------------- | --------------------- |
| Bracken    | short              | Nucleotide             | k-mer             |                       |
| Centrifuge | short              | Nucleotide             | FM-index          |                       |
| DIAMOND    | short/long         | Amino acid             | alignment based   |                       |
| ganon      | short              | Nucleotide             | k-mer             |                       |
| Kaiju      | short              | Amino acid             | alignment based   |                       |
| KMCP       | short              | Nucleotide             | k-mer             |                       |
| Kraken2    | short              | Nucleotide             | k-mer             |                       |
| KrakenUniq | short              | Nucleotide             | k-mer             |                       |
| MALT       | short              | Nucleotide/ Amino acid | alignment based   |                       |
| melon      | long               | Nucleotide             | alignment based   |                       |
| MetaPhlAn  | short              | Nucleotide             | alignment based   |                       |
| MetaCache  | long               | Nucleotide             | k-mer             |                       |
| mOTUs      | short/long         | Nucleotide             | alignment based   |                       |
| sylph      | short/long         | Nucleotide             | k-mer             |                       |

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
