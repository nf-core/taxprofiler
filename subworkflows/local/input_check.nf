//
// Check input samplesheet and get read channels
//

include { EIDO_VALIDATE } from '../../modules/nf-core/modules/eido/validate/main'
include { EIDO_CONVERT } from '../../modules/nf-core/modules/eido/convert/main'

workflow INPUT_CHECK {
    take:
    samplesheet_or_pep_config // file: /path/to/samplesheet.csv or /path/to/pep/config.yaml
    pep_input_base_dir

    main:
    ch_versions = Channel.empty()

    EIDO_VALIDATE ( samplesheet_or_pep_config, file("$projectDir/assets/samplesheet_schema.yaml"), pep_input_base_dir )
    ch_versions = ch_versions.mix(EIDO_VALIDATE.out.versions)

    EIDO_CONVERT ( samplesheet_or_pep_config, "csv", pep_input_base_dir )
    ch_versions = ch_versions.mix(EIDO_CONVERT.out.versions)

    ch_parsed_samplesheet = EIDO_CONVERT.out.samplesheet_converted
        .splitCsv ( header:true, sep:',' )
        .map{

            // Checks not supported by EIDO(?)
            if ( ( it['fastq_1'] != "" || it['fastq_2'] != "" ) && it['fasta'] != "" ) { exit 1, "[nf-core/taxprofiler] ERROR: FastQ and FastA files cannot be specified together in the same library. Check input samplesheet! Check sample: ${it['sample']}" }
            if ( it['fastq_1'] == "" && it['fastq_2'] != "" )                          { exit 1, "[nf-core/taxprofiler] ERROR: Input samplesheet has a missing fastq_1 when fastq_2 is specified. Check sample: ${it['sample']}" }

            single_end = it['fastq_2'] == "" ? true : false
            it['single_end'] = single_end

            [ it ]
        }
        .flatten()
        .branch {
            fasta: it['fasta'] != ''
            nanopore: it['instrument_platform'] == 'OXFORD_NANOPORE'
            fastq: true
        }

    ch_parsed_samplesheet.fastq
        .map { create_fastq_channel(it) }
        .set { fastq }

    ch_parsed_samplesheet.nanopore
        .map { create_fastq_channel(it) }
        .set { nanopore }

    ch_parsed_samplesheet.fasta
        .map { create_fasta_channel(it) }
        .set { fasta }

    emit:
    fastq = fastq ?: []                       // channel: [ val(meta), [ reads ] ]
    nanopore = nanopore ?: []                 // channel: [ val(meta), [ reads ] ]
    fasta = fasta ?: []                       // channel: [ val(meta), fasta ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                     = row.sample
    meta.run_accession          = row.run_accession
    meta.instrument_platform    = row.instrument_platform
    meta.single_end             = row.single_end.toBoolean()
    meta.is_fasta               = false

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (meta.instrument_platform == 'OXFORD_NANOPORE') {
            if (row.fastq_2 != '') {
                exit 1, "ERROR: Please check input samplesheet -> For Oxford Nanopore reads Read 2 FastQ should be empty!\n${row.fastq_2}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1) ] ]
        } else {
            if (!file(row.fastq_2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        }

    }
    return fastq_meta
}// Function to get list of [ meta, fasta ]
def create_fasta_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id                     = row.sample
    meta.run_accession          = row.run_accession
    meta.instrument_platform    = row.instrument_platform
    meta.single_end             = true
    meta.is_fasta               = true

    def array = []
    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> FastA file does not exist!\n${row.fasta}"
    }
    array = [ meta, [ file(row.fasta) ] ]

    return array
}
