//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    parsed_samplesheet = SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .branch { row ->
            fasta: row.fasta != ''
            nanopore: row.instrument_platform == 'OXFORD_NANOPORE'
            fastq: true
        }

    fastq = parsed_samplesheet.fastq
        .map { create_fastq_channel(it) }

    nanopore = parsed_samplesheet.nanopore
        .map { create_fastq_channel(it) }

    fasta = parsed_samplesheet.fasta
        .map { create_fasta_channel(it) }

    emit:
    fastq = fastq ?: []                       // channel: [ val(meta), [ reads ] ]
    nanopore = nanopore ?: []                 // channel: [ val(meta), [ reads ] ]
    fasta = fasta ?: []                       // channel: [ val(meta), fasta ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = row.subMap(['sample', 'run_accession', 'instrument_platform'])
    meta.id         = meta.sample
    meta.single_end = row.single_end.toBoolean()
    meta.is_fasta   = false

    // add path(s) of the fastq file(s) to the meta map
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (meta.single_end) {
        return [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (meta.instrument_platform == 'OXFORD_NANOPORE') {
            if (row.fastq_2 != '') {
                exit 1, "ERROR: Please check input samplesheet -> For Oxford Nanopore reads Read 2 FastQ should be empty!\n${row.fastq_2}"
            }
            return [ meta, [ file(row.fastq_1) ] ]
        } else {
            if (!file(row.fastq_2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
            }
            return [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        }
    }
}

// Function to get list of [ meta, fasta ]
def create_fasta_channel(LinkedHashMap row) {
    def meta        = row.subMap(['sample', 'run_accession', 'instrument_platform'])
    meta.id         = meta.sample
    meta.single_end = true
    meta.is_fasta   = true

    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> FastA file does not exist!\n${row.fasta}"
    }
    return [ meta, [ file(row.fasta) ] ]
}
