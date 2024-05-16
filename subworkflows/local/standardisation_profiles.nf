//
// Standardise output files e.g. aggregation
//

include { TAXPASTA_MERGE                                                        } from '../../modules/nf-core/taxpasta/merge/main'
include { TAXPASTA_STANDARDISE                                                  } from '../../modules/nf-core/taxpasta/standardise/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS                                         } from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { KAIJU_KAIJU2TABLE as KAIJU_KAIJU2TABLE_COMBINED                       } from '../../modules/nf-core/kaiju/kaiju2table/main'
include { KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_KRAKEN     } from '../../modules/nf-core/krakentools/combinekreports/main'
include { KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE } from '../../modules/nf-core/krakentools/combinekreports/main'
include { METAPHLAN_MERGEMETAPHLANTABLES                                        } from '../../modules/nf-core/metaphlan/mergemetaphlantables/main'
include { MOTUS_MERGE                                                           } from '../../modules/nf-core/motus/merge/main'
include { GANON_TABLE                                                           } from '../../modules/nf-core/ganon/table/main'

// Custom Functions

/**
* Group all profiles per reference database.
*
* @param ch_profiles A channel containing pairs of a meta map and the report of
*   a given profiler, where meta must contain a key `db_name`.
* @return A channel with one element per reference database. Each element is a
*   pair of a meta map with an `id` key and all corresponding profiles.
*/
def groupProfiles(ch_profiles, groupTupleOptions = [:]) {
    return ch_profiles
        .map { meta, profile -> [meta.db_name, profile] }
        .groupTuple(groupTupleOptions)
        .map { db_name, profiles -> [[id: db_name], profiles] }
}

/**
* Combine profiles with their corresponding reference database, then separate into two channels.
*
* The combined results are returned on multiple channels, where the element
* position for the profiles in one channel is the same as the position of the
* corresponding database element in the other channel.
*
* @param ch_profiles A channel containing pairs of a meta map with an `id` key
*   for a reference database, and all the corresponding profiling reports.
* @param ch_database A channel containing pairs of a database meta map and the
*   database itself.
* @return A multiMap'ed output channel with two sub channels, one with the
*   profiles (`profile`) and the other with the corresponding database (`db`).
*/
def combineProfilesWithDatabase(ch_profiles, ch_database) {
    return ch_profiles
        .map { meta, profile -> [meta.id, meta, profile] }
        .combine(ch_database.map { db_meta, db -> [db_meta.db_name, db] }, by: 0)
        .multiMap {
            key, meta, profile, db ->
                profile: [meta, profile]
                db: db
        }
}

workflow STANDARDISATION_PROFILES {
    take:
    classifications
    profiles
    databases
    motu_version

    main:
    ch_versions            = Channel.empty()
    ch_multiqc_files       = Channel.empty()

    //Taxpasta standardisation
    ch_prepare_for_taxpasta = profiles
                            .map {
                                    meta, profile ->
                                        def meta_new = [:]
                                        meta_new.tool = meta.tool == 'malt' ? 'megan6' : meta.tool
                                        meta_new.db_name = meta.db_name
                                        [meta_new, profile]
                            }
                            .groupTuple ()
                            .map {
                                meta, profiles ->
                                    meta = meta + [
                                        tool: meta.tool == 'kraken2-bracken' ? 'kraken2' : meta.tool, // replace to get the right output-format description
                                        id: meta.tool == 'kraken2-bracken' ? "${meta.db_name}-bracken" : "${meta.db_name}" // append so to disambiguate when we have same databases for kraken2 step of bracken, with normal bracken
                                    ]
                                [ meta, profiles.flatten() ]
                            }

    ch_taxpasta_tax_dir = params.taxpasta_taxonomy_dir ? Channel.fromPath(params.taxpasta_taxonomy_dir, checkIfExists: true).collect() : []

    ch_input_for_taxpasta = ch_prepare_for_taxpasta
                        .branch {
                            meta, profile ->
                                merge:      profile.size() > 1
                                standardise: true
                        }


    TAXPASTA_MERGE       (ch_input_for_taxpasta.merge      , ch_taxpasta_tax_dir, [])
    ch_versions = ch_versions.mix( TAXPASTA_MERGE.out.versions.first() )
    TAXPASTA_STANDARDISE (ch_input_for_taxpasta.standardise, ch_taxpasta_tax_dir    )
    ch_version = ch_versions.mix( TAXPASTA_STANDARDISE.out.versions.first() )



    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles
        .branch {
            bracken: it[0]['tool'] == 'bracken'
            centrifuge: it[0]['tool'] == 'centrifuge'
            ganon: it[0]['tool'] == 'ganon'
            kmcp: it [0]['tool'] == 'kmcp'
            kraken2: it[0]['tool'] == 'kraken2' || it[0]['tool'] == 'kraken2-bracken'
            metaphlan: it[0]['tool'] == 'metaphlan'
            motus: it[0]['tool'] == 'motus'
            unknown: true
        }

    ch_input_classifications = classifications
        .branch {
            kaiju: it[0]['tool'] == 'kaiju'
            unknown: true
        }

    ch_input_databases = databases
        .branch {
            motus: it[0]['tool'] == 'motus'
            kaiju: it[0]['tool'] == 'kaiju'
            unknown: true
        }

    /*
        Standardise and aggregate
    */

    // Bracken

    ch_profiles_for_bracken = groupProfiles(ch_input_profiles.bracken)

    BRACKEN_COMBINEBRACKENOUTPUTS ( ch_profiles_for_bracken )

    // CENTRIFUGE

    // Collect and replace id for db_name for prefix
    // Have to sort by size to ensure first file actually has hits otherwise
    // the script fails
    ch_profiles_for_centrifuge = groupProfiles(
        ch_input_profiles.centrifuge,
        [sort: { -it.size() }]
    )

    KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE ( ch_profiles_for_centrifuge )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE.out.versions )

    // Kaiju

    // Collect and replace id for db_name for prefix
    ch_profiles_for_kaiju = groupProfiles(ch_input_classifications.kaiju)

    ch_input_for_kaiju2tablecombine = combineProfilesWithDatabase(ch_profiles_for_kaiju, ch_input_databases.kaiju)

    KAIJU_KAIJU2TABLE_COMBINED ( ch_input_for_kaiju2tablecombine.profile, ch_input_for_kaiju2tablecombine.db, params.kaiju_taxon_rank)
    ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE_COMBINED.out.summary )
    ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE_COMBINED.out.versions )

    // Kraken2

    // Collect and replace id for db_name for prefix
    // Have to sort by size to ensure first file actually has hits otherwise
    // the script fails
    ch_profiles_for_kraken2 = groupProfiles(
        ch_input_profiles.kraken2
        .map { meta, profile ->
            // Replace database name, to get the right output description.
            def db_name = meta.tool == 'kraken2-bracken' ? "${meta.db_name}-bracken" : "${meta.db_name}"
            return [meta + [db_name: db_name], profile]
        },
        [sort: { -it.size() }]
    )

    KRAKENTOOLS_COMBINEKREPORTS_KRAKEN ( ch_profiles_for_kraken2 )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKENTOOLS_COMBINEKREPORTS_KRAKEN.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_COMBINEKREPORTS_KRAKEN.out.versions )

    // MetaPhlAn

    ch_profiles_for_metaphlan = groupProfiles(ch_input_profiles.metaphlan)

    METAPHLAN_MERGEMETAPHLANTABLES ( ch_profiles_for_metaphlan )
    ch_multiqc_files = ch_multiqc_files.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.txt )
    ch_versions = ch_versions.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.versions )

    // mOTUs

    // mOTUs has a 'single' database, and cannot create custom ones.
    // Therefore removing db info here, and publish merged at root mOTUs results
    // directory

    ch_profiles_for_motus = groupProfiles(ch_input_profiles.motus)

    ch_input_for_motusmerge = combineProfilesWithDatabase(ch_profiles_for_motus, ch_input_databases.motus)

    MOTUS_MERGE ( ch_input_for_motusmerge.profile, ch_input_for_motusmerge.db, motu_version )
    ch_versions = ch_versions.mix( MOTUS_MERGE.out.versions )

    // Ganon

    ch_profiles_for_ganon = groupProfiles(ch_input_profiles.ganon)

    GANON_TABLE ( ch_profiles_for_ganon )
    ch_multiqc_files = ch_multiqc_files.mix( GANON_TABLE.out.txt )
    ch_versions = ch_versions.mix( GANON_TABLE.out.versions )

    emit:
    taxpasta = TAXPASTA_MERGE.out.merged_profiles
    versions = ch_versions
    mqc      = ch_multiqc_files
}
