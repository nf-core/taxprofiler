//
// Create Krona visualizations
//

include { KRAKENTOOLS_KREPORT2KRONA } from '../../modules/nf-core/modules/krakentools/kreport2krona/main'
include { KRONA_CLEANUP             } from '../../modules/local/krona_cleanup'
include { KRONA_KTIMPORTTEXT        } from '../../modules/nf-core/modules/krona/ktimporttext/main'

workflow VISUALIZATION_KRONA {
    take:
    profiles

    main:
    ch_krona_text = Channel.empty()
    ch_krona_html = Channel.empty()
    ch_versions = Channel.empty()

    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles
        .branch {
            kraken2: it[0]['tool'] == 'kraken2'
            unknown: true
        }

    /*
        Convert Kraken2 formatted reports into Krona text files
    */
    ch_kraken_reports = ch_input_profiles.kraken2
    KRAKENTOOLS_KREPORT2KRONA ( ch_kraken_reports )
    ch_krona_text = ch_krona_text.mix( KRAKENTOOLS_KREPORT2KRONA.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_KREPORT2KRONA.out.versions.first() )

    /*
        Remove taxonomy level annotations from the Krona text files
    */
    KRONA_CLEANUP( ch_krona_text )
    ch_cleaned_krona_text = KRONA_CLEANUP.out.txt
    ch_versions = ch_versions.mix( KRONA_CLEANUP.out.versions.first() )

    /*
        Convert Krona text files into html Krona visualizations
    */
    ch_krona_text_for_import = ch_cleaned_krona_text
        .map{[[id: it[0]['db_name']], it[1]]}
        .groupTuple()
    KRONA_KTIMPORTTEXT( ch_krona_text_for_import )
    ch_krona_html = ch_krona_html.mix( KRONA_KTIMPORTTEXT.out.html )
    ch_versions = ch_versions.mix( KRONA_KTIMPORTTEXT.out.versions.first() )

    emit:
    html = ch_krona_html
    versions = ch_versions
}
