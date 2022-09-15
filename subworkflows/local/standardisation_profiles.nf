//
// Standardise output files e.g. aggregation
//

include { KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_CENTRIFUGE } from '../../modules/nf-core/modules/krakentools/combinekreports/main'
include { KAIJU_KAIJU2TABLE                                                     } from '../../modules/nf-core/modules/kaiju/kaiju2table/main'
include { KRAKENTOOLS_COMBINEKREPORTS                                           } from '../../modules/nf-core/modules/krakentools/combinekreports/main'
include { MOTUS_MERGE                                                           } from '../../modules/nf-core/modules/motus/merge/main'

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
            centrifuge: it[0]['tool'] == 'centrifuge'
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

    // mOTUs has a 'single' database, and cannot create custom ones.
    // Therefore removing db info here, and publish merged at root mOTUs results
    // directory
    MOTUS_MERGE ( ch_input_profiles.motus.map{it[1]}.collect(), ch_input_databases.motus.map{it[1]}, motu_version, params.generate_biom_output )
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
