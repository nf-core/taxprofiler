//
// Create Krona visualizations
//

include { MOTUS_MERGE } from '../../modules/nf-core/modules/motus/merge/main'

workflow STANDARDISATION_PROFILES {
    take:
    classifications
    profiles
    databases
    motu_version

    main:
    ch_standardised_tables = Channel.empty()
    ch_versions            = Channel.empty()

    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles
        .branch {
            motus: it[0]['tool'] == 'motus'
            unknown: true
        }

    ch_input_classifications = classifications
        .branch {
            unknown: true
        }

    ch_input_databases = databases
        .branch {
            motus: it[0]['tool'] == 'motus'
            unknown: true
        }

    /*
        Standardise and aggregate
    */

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
    tables = ch_standardised_tables
    versions = ch_versions
}
