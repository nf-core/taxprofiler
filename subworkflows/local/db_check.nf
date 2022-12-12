//
// Check input samplesheet and get read channels
//

include { UNTAR            } from '../../modules/nf-core/untar/main'

workflow DB_CHECK {
    take:
    dbsheet // file: /path/to/dbsheet.csv

    main:
    ch_versions = Channel.empty()

    // TODO: make database sheet check
    // Checks:
    // 1) no duplicates,

    parsed_samplesheet = Channel.fromPath(dbsheet)
        .splitCsv ( header:true, sep:',' )
        .map {
            validate_db_sheet(it)
            create_db_channels(it)
        }

    ch_dbs_for_untar = parsed_samplesheet
        .branch {
            untar: it[1].toString().endsWith(".tar.gz")
            skip: true
        }

    // TODO Filter to only run UNTAR on DBs of tools actually using?
    // TODO make optional whether to save
    UNTAR ( ch_dbs_for_untar.untar )
    ch_versions = ch_versions.mix(UNTAR.out.versions.first())

    ch_final_dbs = ch_dbs_for_untar.skip.mix( UNTAR.out.untar )

    emit:
    dbs = ch_final_dbs                        // channel: [ val(meta), [ db ] ]
    versions = ch_versions                   // channel: [ versions.yml ]
}

def validate_db_sheet(LinkedHashMap row){

        // check minimum number of columns
        if (row.size() < 4) exit 1, "[nf-core/taxprofiler] error:  Invalid database input sheet - malformed row (e.g. missing column). See documentation for more information. Error in: ${row}, "

        // all columns there
        def expected_headers = ['tool', 'db_name', 'db_params', 'db_path']
        if ( !row.keySet().containsAll(expected_headers) ) exit 1, "[nf-core/taxprofiler] error: Invalid database input sheet - malformed column names. Please check input TSV. Column names should be: ${expected_keys.join(", ")}"

        // valid tools specified// TIFNISIH LIST
        def expected_tools = [ "bracken", "centrifuge", "diamond", "kaiju",  "kraken2", "malt", "metaphlan3"  ]

        // detect quotes in params
        if ( row.db_params.contains('"') ) exit 1, "[nf-core/taxprofiler] error: Invalid database db_params entry. No quotes allowed. Error in: ${row}"
        if ( row.db_params.contains("'") ) exit 1, "[nf-core/taxprofiler] error: Invalid database db_params entry. No quotes allowed. Error in: ${row}"

}

def create_db_channels(LinkedHashMap row) {
    def meta = [:]
    meta.tool             = row.tool
    meta.db_name          = row.db_name
    meta.db_params        = row.db_params

    def array = []
    if (!file(row.db_path, type: 'dir').exists()) {
        exit 1, "ERROR: Please check input samplesheet -> database could not be found!\n${row.db_path}"
    }
    array = [ meta, file(row.db_path) ]

    return array
}


