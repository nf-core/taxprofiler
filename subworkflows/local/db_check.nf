//
// Check input samplesheet and get read channels
//

include { DATABASE_CHECK } from '../../modules/local/database_check'

workflow DB_CHECK {
    take:
    dbsheet // file: /path/to/dbsheet.csv

    main:

    // TODO: make database sheet check
    parsed_samplesheet = DATABASE_CHECK ( dbsheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .dump(tag: "db_split_csv_out")
        .map { create_db_channels(it) }
        .dump(tag: "db_channel_prepped")
        .set{ dbs }

    emit:
    dbs                                       // channel: [ val(meta), [ db ] ]
    versions = DATABASE_CHECK.out.versions // channel: [ versions.yml ]
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
