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
include { KMCP_PROFILE                                                          } from '../../modules/nf-core/kmcp/profile/main'
include { GANON_TABLE                                                           } from '../../modules/nf-core/ganon/table/main'

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
                            .filter {
                                meta, report ->
                                // TODO: add tool to taxpasta
                                    if ( meta['tool'] == 'kmcp' ) log.warn "[nf-core/taxprofiler] kmcp is not yet supported in Taxpasta. Skipping kmcp profile for sample ${meta.id}."
                                    meta['tool'] != 'kmcp'
                             }
                            .map {
                                    meta, profile ->
                                        def meta_new = [:]
                                        meta_new.id = meta.db_name
                                        meta_new.tool = meta.tool == 'metaphlan' ? 'metaphlan' : meta.tool == 'malt' ? 'megan6' : meta.tool
                                        [meta_new, profile]
                            }
                            .groupTuple ()
                            .map { [ it[0], it[1].flatten() ] }

    ch_taxpasta_tax_dir = params.taxpasta_taxonomy_dir ? Channel.fromPath(params.taxpasta_taxonomy_dir, checkIfExists: true).collect() : []

    ch_input_for_taxpasta = ch_prepare_for_taxpasta
                        .branch {
                            meta, profile ->
                                merge:      profile.size() > 1
                                standardise: true
                        }


    TAXPASTA_MERGE       (ch_input_for_taxpasta.merge      , ch_taxpasta_tax_dir, [])
    TAXPASTA_STANDARDISE (ch_input_for_taxpasta.standardise, ch_taxpasta_tax_dir    )

    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles
        .branch {
            bracken: it[0]['tool'] == 'bracken'
            centrifuge: it[0]['tool'] == 'centrifuge'
            kraken2: it[0]['tool'] == 'kraken2'
            metaphlan: it[0]['tool'] == 'metaphlan'
            motus: it[0]['tool'] == 'motus'
            kmcp: it [0]['tool'] == 'kmcp'
            ganon: it[0]['tool'] == 'ganon'
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

    ch_profiles_for_bracken = ch_input_profiles.bracken
                            .map { [it[0]['db_name'], it[1]] }
                            .groupTuple()
                            .map {
                                [[id:it[0]], it[1]]
                            }

    BRACKEN_COMBINEBRACKENOUTPUTS ( ch_profiles_for_bracken )

    // CENTRIFUGE

    // Collect and replace id for db_name for prefix
    // Have to sort by size to ensure first file actually has hits otherwise
    // the script fails
    ch_profiles_for_centrifuge = ch_input_profiles.centrifuge
                                .map { [it[0]['db_name'], it[1]] }
                                .groupTuple(sort: {-it.size()} )
                                .map {
                                    [[id:it[0]], it[1]]
                                }


    KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE ( ch_profiles_for_centrifuge )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE.out.versions )

    // Kaiju

    // Collect and replace id for db_name for prefix
    ch_profiles_for_kaiju = ch_input_classifications.kaiju
                                .map { [it[0]['db_name'], it[1]] }
                                .groupTuple()
                                .map {
                                    [[id:it[0]], it[1]]
                                }

    ch_input_for_kaiju2tablecombine = ch_profiles_for_kaiju
                                        .map { meta, profile -> [meta.id, meta, profile] }
                                        .combine(ch_input_databases.kaiju.map{meta, db -> [meta.db_name, meta, db]}, by: 0)
                                        .multiMap {
                                            key, meta, profile, db_meta, db ->
                                                profile: [meta, profile]
                                                db: db
                                        }

    KAIJU_KAIJU2TABLE_COMBINED ( ch_input_for_kaiju2tablecombine.profile, ch_input_for_kaiju2tablecombine.db, params.kaiju_taxon_rank)
    ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE_COMBINED.out.summary )
    ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE_COMBINED.out.versions )

    // Kraken2

    // Collect and replace id for db_name for prefix
    // Have to sort by size to ensure first file actually has hits otherwise
    // the script fails
    ch_profiles_for_kraken2 = ch_input_profiles.kraken2
                                .map { [it[0]['db_name'], it[1]] }
                                .groupTuple(sort: {-it.size()} )
                                .map {
                                    [[id:it[0]], it[1]]
                                }

    KRAKENTOOLS_COMBINEKREPORTS_KRAKEN ( ch_profiles_for_kraken2 )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKENTOOLS_COMBINEKREPORTS_KRAKEN.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_COMBINEKREPORTS_KRAKEN.out.versions )

    // MetaPhlAn

    ch_profiles_for_metaphlan = ch_input_profiles.metaphlan
                            .map { [it[0]['db_name'], it[1]] }
                            .groupTuple()
                            .map {
                                [[id:it[0]], it[1]]
                            }

    METAPHLAN_MERGEMETAPHLANTABLES ( ch_profiles_for_metaphlan )
    ch_multiqc_files = ch_multiqc_files.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.txt )
    ch_versions = ch_versions.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.versions )

    // mOTUs

    // mOTUs has a 'single' database, and cannot create custom ones.
    // Therefore removing db info here, and publish merged at root mOTUs results
    // directory

    ch_profiles_for_motus = ch_input_profiles.motus
                                .map { [it[0]['db_name'], it[1]] }
                                .groupTuple()
                                .map {
                                    [[id:it[0]], it[1]]
                                }

    ch_input_for_motusmerge = ch_profiles_for_motus
                                        .map { meta, profile -> [meta.id, meta, profile] }
                                        .combine(ch_input_databases.motus.map{meta, db -> [meta.db_name, meta, db]}, by: 0)
                                        .multiMap {
                                            key, meta, profile, db_meta, db ->
                                                profile: [meta, profile]
                                                db: db
                                        }

    MOTUS_MERGE ( ch_input_for_motusmerge.profile, ch_input_for_motusmerge.db, motu_version )
    ch_versions = ch_versions.mix( MOTUS_MERGE.out.versions )

    // Ganon

    ch_profiles_for_ganon = ch_input_profiles.ganon
                            .map { [it[0]['db_name'], it[1]] }
                            .groupTuple()
                            .map {
                                [[id:it[0]], it[1]]
                            }

    GANON_TABLE ( ch_profiles_for_ganon )
    ch_multiqc_files = ch_multiqc_files.mix( GANON_TABLE.out.txt )
    ch_versions = ch_versions.mix( GANON_TABLE.out.versions )

    emit:
    taxpasta = TAXPASTA_MERGE.out.merged_profiles
    versions = ch_versions
    mqc      = ch_multiqc_files
}
