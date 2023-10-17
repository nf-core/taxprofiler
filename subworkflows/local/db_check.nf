//
// Check input samplesheet and get read channels
//

include { UNTAR } from '../../modules/nf-core/untar/main'

workflow DB_CHECK {
    take:
    dbsheet // file: /path/to/dbsheet.csv

    main:
    ch_versions = Channel.empty()
    ch_dbs_for_untar = Channel.empty()
    ch_final_dbs = Channel.empty()

    /*
        Validation
    */

    // Special check to header exists
    Channel.fromPath(dbsheet)
            .splitCsv ( header:false, sep:',' )
            .first()
            .map {
                if ( it[0] != 'tool' && it[1] != 'db_name' && it[2] != 'db_params' && it[2] != 'db_path' ) error ("[nf-core/taxprofiler] ERROR: database sheet is missing header. Please see nf-core/taxprofiler documentation for database sheet specifications.")
            }

    // Normal checks for within-row validity, so can be moved to separate functions
    parsed_samplesheet = Channel.fromPath(dbsheet)
        .splitCsv ( header:true, sep:',' )
        .map { row ->
            validate_db_rows(row)
            return [ row.subMap(['tool', 'db_name', 'db_params']), file(row.db_path) ]
        }

    // Special check to check _between_ rows, for which we must group rows together
    // Note: this will run in parallel to within-row validity, but we can assume this will run faster thus will fail first
    Channel.fromPath(dbsheet)
            .splitCsv ( header:true, sep:',' )
            .map {[it.tool, it.db_name] }
            .groupTuple()
            .map {
                tool, db_name ->
                    def unique_names = db_name.unique(false)
                    if ( unique_names.size() < db_name.size() ) error("[nf-core/taxprofiler] ERROR: Each database for a tool must have a unique name, duplicated detected. Tool: ${tool}, Database name: ${unique_names}")
            }

    /*
        Preparation
    */

    // Decompress
    ch_dbs_for_untar = parsed_samplesheet
        .branch { db_meta, db ->
            untar: db.name.endsWith(".tar.gz")
            skip: true
        }

    // Filter the channel to untar only those databases for tools that are selected to be run by the user.
    ch_input_untar = ch_dbs_for_untar.untar
        .filter { db_meta, db -> params["run_${db_meta.tool}"] }

    UNTAR (ch_input_untar)
    ch_versions = ch_versions.mix(UNTAR.out.versions.first())
    ch_final_dbs = ch_dbs_for_untar.skip.mix( UNTAR.out.untar )

    emit:
    dbs = ch_final_dbs                        // channel: [ val(meta), [ db ] ]
    versions = ch_versions                   // channel: [ versions.yml ]
}

def validate_db_rows(LinkedHashMap row) {

    // check minimum number of columns
    if (row.size() < 4) error("[nf-core/taxprofiler] ERROR: Invalid database input sheet - malformed row (e.g. missing column). See documentation for more information. Error in: ${row}")

    // all columns there
    def expected_headers = ['tool', 'db_name', 'db_params', 'db_path']
    if ( !row.keySet().containsAll(expected_headers) ) error("[nf-core/taxprofiler] ERROR: Invalid database input sheet - malformed column names. Please check input TSV. Column names should be: ${expected_headers.join(", ")}")

    // valid tools specified
    def expected_tools = [ "bracken", "centrifuge", "diamond", "ganon", "kaiju", "kmcp", "kraken2", "krakenuniq", "malt", "metaphlan", "motus", "sourmash" ]

    if ( !expected_tools.contains(row.tool) ) error("[nf-core/taxprofiler] ERROR: Invalid tool name. Please see documentation for all supported profilers. Error in: ${row}")

    // detect quotes in params
    if ( row.db_params.contains('"') ) error("[nf-core/taxprofiler] ERROR: Invalid database db_params entry. No quotes allowed. Error in: ${row}")
    if ( row.db_params.contains("'") ) error("[nf-core/taxprofiler] ERROR: Invalid database db_params entry. No quotes allowed. Error in: ${row}")

    // check if any form of bracken params, that it must have `;`
    if ( row.tool == 'bracken' && row.db_params && !row.db_params.contains(";") )  error("[nf-core/taxprofiler] ERROR: Invalid database db_params entry. Bracken requires a semi-colon if passing parameter. Error in: ${row}")
    if ( row.tool == 'kmcp' && row.db_params && row.db_params.matches(".*;\$"))  error ("[nf-core/taxprofiler] ERROR: Invalid database db_params entry. KMCP only requires a semi-colon if passing arguments to KMCP profile, in cases of which the arguments should go after the semi-colon. Error in: ${row}")

    // ensure that the database directory exists
    if (!file(row.db_path, type: 'dir').exists()) error("ERROR: Please check input samplesheet -> database path could not be found!\n${row.db_path}")

}
