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
        .dump(tag: "input_split_csv_out")
        .branch {
            fasta: it['fasta'] != ''
            fastq: true
        }

    parsed_samplesheet.fastq
        .map { create_fastq_channels(it) }
        .dump(tag: "fastq_channel_init")
        .set { fastq }

    parsed_samplesheet.fasta
        .map { create_fasta_channels(it) }
        .dump(tag: "fasta_channel_init")
        .set { fasta }

    emit:
    fastq                                     // channel: [ val(meta), [ reads ] ]
    fasta                                     // channel: [ val(meta), fasta ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                     = row.sample
    meta.run_accession          = row.run_accession
    meta.instrument_platform    = row.instrument_platform
    meta.single_end             = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

// Function to get list of [ meta, fasta ]
def create_fasta_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id                     = row.sample
    meta.run_accession          = row.run_accession
    meta.instrument_platform    = row.instrument_platform
    meta.single_end             = true

    def array = []
    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> FastA file does not exist!\n${row.fasta}"
    }
    array = [ meta, [ file(row.fasta) ] ]

    return array
}
