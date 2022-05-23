//
// Check input samplesheet and get read channels
//

include { DATABASE_CHECK   } from '../../modules/local/database_check'
include { UNTAR            } from '../../modules/nf-core/modules/untar/main'
include { MOTUS_DOWNLOADDB } from '../../modules/nf-core/modules/motus/downloaddb/main'

workflow DB_CHECK {
    take:
    dbsheet // file: /path/to/dbsheet.csv

    main:

    // TODO: make database sheet check
    // Checks:
    // 1) no duplicates,
    // 2) args do not have quotes, e.g. just `,,` and NOT `,"",`
    parsed_samplesheet = DATABASE_CHECK ( dbsheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_db_channels(it) }

    // Download database for mOTUs
    if( params.run_motus ){
        check_motus_db =
            parsed_samplesheet.filter{ it[0].tool == "motus" }
                .ifEmpty{[]}
        if( params.download_motus_db ){
            MOTUS_DOWNLOADDB( params.motus_downloaddb_script )
            check_motus_db = MOTUS_DOWNLOADDB.out.db
                .map{[
                    [tool: "motus", db_name: "db_mOTU", db_params: ''],
                    it
                ]}
        }
        parsed_samplesheet = parsed_samplesheet.mix(check_motus_db)
    }

    ch_dbs_for_untar = parsed_samplesheet
        .branch {
            untar: it[1].toString().endsWith(".tar.gz")
            skip: true
        }

    // TODO Filter to only run UNTAR on DBs of tools actually using?
    // TODO make optional whether to save
    UNTAR ( ch_dbs_for_untar.untar )

    ch_final_dbs = ch_dbs_for_untar.skip.mix( UNTAR.out.untar )

    emit:
    dbs = ch_final_dbs                        // channel: [ val(meta), [ db ] ]
    versions = DATABASE_CHECK.out.versions.mix(UNTAR.out.versions.first()) // channel: [ versions.yml ]
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
