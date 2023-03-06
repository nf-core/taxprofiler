//
// Standardise output files e.g. aggregation
//

include { BRACKEN_COMBINEBRACKENOUTPUTS                                         } from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { KAIJU_KAIJU2TABLE as KAIJU_KAIJU2TABLE_COMBINED                       } from '../../modules/nf-core/kaiju/kaiju2table/main'
include { KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_KRAKEN     } from '../../modules/nf-core/krakentools/combinekreports/main'
include { KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE } from '../../modules/nf-core/krakentools/combinekreports/main'
include { METAPHLAN3_MERGEMETAPHLANTABLES                                       } from '../../modules/nf-core/metaphlan3/mergemetaphlantables/main'
include { MOTUS_MERGE                                                           } from '../../modules/nf-core/motus/merge/main'
include { TAXPASTA_MERGE                                                        } from '../../modules/nf-core/taxpasta/merge/main'

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
    ch_input_for_taxpasta = profiles
                            .map {
                                    meta, profile ->
                                        def meta_new = [:]
                                        meta_new.id = meta.db_name
                                        meta_new.tool = meta.tool == 'metaphlan3' ? 'metaphlan' : meta.tool == 'malt' ? 'megan6' : meta.tool
                                        [meta_new, profile]
                            }
                            .groupTuple ()
                            .map { [ it[0], it[1].flatten() ] }

    TAXPASTA_MERGE (ch_input_for_taxpasta, [], [])

    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles
        .branch {
            bracken: it[0]['tool'] == 'bracken'
            centrifuge: it[0]['tool'] == 'centrifuge'
            kraken2: it[0]['tool'] == 'kraken2'
            metaphlan3: it[0]['tool'] == 'metaphlan3'
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

    KAIJU_KAIJU2TABLE_COMBINED ( ch_profiles_for_kaiju, ch_input_databases.kaiju.map{it[1]}, params.kaiju_taxon_rank)
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

    // MetaPhlAn3

    ch_profiles_for_metaphlan3 = ch_input_profiles.metaphlan3
                            .map { [it[0]['db_name'], it[1]] }
                            .groupTuple()
                            .map {
                                [[id:it[0]], it[1]]
                            }

    METAPHLAN3_MERGEMETAPHLANTABLES ( ch_profiles_for_metaphlan3 )
    ch_multiqc_files = ch_multiqc_files.mix( METAPHLAN3_MERGEMETAPHLANTABLES.out.txt )
    ch_versions = ch_versions.mix( METAPHLAN3_MERGEMETAPHLANTABLES.out.versions )

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

    MOTUS_MERGE ( ch_profiles_for_motus, ch_input_databases.motus.map{it[1]}, motu_version )
    ch_versions = ch_versions.mix( MOTUS_MERGE.out.versions )

    emit:
    taxpasta = TAXPASTA_MERGE.out.merged_profiles
    versions = ch_versions
    mqc      = ch_multiqc_files
}
