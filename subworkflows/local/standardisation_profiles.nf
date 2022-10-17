//
// Standardise output files e.g. aggregation
//

include { KAIJU_KAIJU2TABLE                                                     } from '../../modules/nf-core/kaiju/kaiju2table/main'
include {
    KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE;
    KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_BRACKEN;
    KRAKENTOOLS_COMBINEKREPORTS
} from '../../modules/nf-core/krakentools/combinekreports/main'
include { METAPHLAN3_MERGEMETAPHLANTABLES                                       } from '../../modules/nf-core/metaphlan3/mergemetaphlantables/main'
include { MOTUS_MERGE                                                           } from '../../modules/nf-core/motus/merge/main'

workflow STANDARDISATION_PROFILES {
    take:
    classifications
    profiles
    databases
    motu_version

    main:
    ch_standardised_tables = Channel.empty()
    ch_versions            = Channel.empty()
    ch_multiqc_files       = Channel.empty()

    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles
        .branch {
            motus: it[0]['tool'] == 'motus'
            kraken2: it[0]['tool'] == 'kraken2'
            bracken: it[0]['tool'] == 'bracken'
            centrifuge: it[0]['tool'] == 'centrifuge'
            metaphlan3: it[0]['tool'] == 'metaphlan3'
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
    ch_standardised_tables = ch_standardised_tables.mix( KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE.out.txt )
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

    KAIJU_KAIJU2TABLE ( ch_profiles_for_kaiju, ch_input_databases.kaiju.map{it[1]}, params.kaiju_taxon_rank)
    ch_standardised_tables = ch_standardised_tables.mix( KAIJU_KAIJU2TABLE.out.summary )
    ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE.out.summary )
    ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE.out.versions )

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

    KRAKENTOOLS_COMBINEKREPORTS ( ch_profiles_for_kraken2 )
    ch_standardised_tables = ch_standardised_tables.mix( KRAKENTOOLS_COMBINEKREPORTS.out.txt )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKENTOOLS_COMBINEKREPORTS.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_COMBINEKREPORTS.out.versions )

    // Bracken

    // Collect and replace id for db_name for prefix
    // Have to sort by size to ensure first file actually has hits otherwise
    // the script fails
    ch_profiles_for_bracken = ch_input_profiles.bracken
                                .map { [it[0]['db_name'], it[1]] }
                                .groupTuple(sort: {-it.size()} )
                                .map {
                                    [[id:it[0]], it[1]]
                                }

    KRAKENTOOLS_COMBINEKREPORTS_BRACKEN ( ch_profiles_for_bracken )
    ch_standardised_tables = ch_standardised_tables.mix( KRAKENTOOLS_COMBINEKREPORTS_BRACKEN.out.txt )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKENTOOLS_COMBINEKREPORTS_BRACKEN.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_COMBINEKREPORTS_BRACKEN.out.versions )

    // MetaPhlAn3

    ch_profiles_for_metaphlan3 = ch_input_profiles.metaphlan3
                            .map { [it[0]['db_name'], it[1]] }
                            .groupTuple()
                            .map {
                                [[id:it[0]], it[1]]
                            }

    METAPHLAN3_MERGEMETAPHLANTABLES ( ch_profiles_for_metaphlan3 )
    ch_standardised_tables = ch_standardised_tables.mix( METAPHLAN3_MERGEMETAPHLANTABLES.out.txt )
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

    if ( params.generate_biom_output ) {
        ch_standardised_tables = ch_standardised_tables.mix ( MOTUS_MERGE.out.biom )
    } else {
        ch_standardised_tables = ch_standardised_tables.mix ( MOTUS_MERGE.out.txt )
    }
    ch_versions = ch_versions.mix( MOTUS_MERGE.out.versions )

    emit:
    tables   = ch_standardised_tables
    versions = ch_versions
    mqc      = ch_multiqc_files
}
