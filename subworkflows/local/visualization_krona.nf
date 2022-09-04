//
// Create Krona visualizations
//

include { MEGAN_RMA2INFO as MEGAN_RMA2INFO_KRONA      } from '../../modules/nf-core/modules/megan/rma2info/main'
include { KAIJU_KAIJU2KRONA         } from '../../modules/nf-core/modules/kaiju/kaiju2krona/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../modules/nf-core/modules/krakentools/kreport2krona/main'
include { KRONA_CLEANUP             } from '../../modules/local/krona_cleanup'
include { KRONA_KTIMPORTTEXT        } from '../../modules/nf-core/modules/krona/ktimporttext/main'
include { KRONA_KTIMPORTTAXONOMY    } from '../../modules/nf-core/modules/krona/ktimporttaxonomy/main'
include { GUNZIP                    } from '../../modules/nf-core/modules/gunzip/main'

workflow VISUALIZATION_KRONA {
    take:
    classifications
    profiles
    databases

    main:
    ch_krona_text = Channel.empty()
    ch_krona_html = Channel.empty()
    ch_versions = Channel.empty()

    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles
        .branch {
            centrifuge: it[0]['tool'] == 'centrifuge'
            kraken2: it[0]['tool'] == 'kraken2'
            unknown: true
        }
    ch_input_classifications = classifications
        .branch {
            kaiju: it[0]['tool'] == 'kaiju'
            malt: it[0]['tool'] == 'malt'
            unknown: true
        }

    /*
        Convert Kraken2 formatted reports into Krona text files
    */
    ch_kraken_reports = ch_input_profiles.kraken2
        .mix( ch_input_profiles.centrifuge )
    KRAKENTOOLS_KREPORT2KRONA ( ch_kraken_reports )
    ch_krona_text = ch_krona_text.mix( KRAKENTOOLS_KREPORT2KRONA.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_KREPORT2KRONA.out.versions.first() )

    /*
        Combine Kaiju profiles with their databases
    */
    ch_input_for_kaiju2krona = ch_input_classifications.kaiju
        .map{ [it[0]['db_name'], it[0], it[1]] }
        .combine( databases.map{ [it[0]['db_name'], it[1]] }, by: 0 )
        .multiMap{
            it ->
                profiles: [it[1], it[2]]
                db: it[3]
        }

    /*
        Convert Kaiju formatted reports into Krona text files
    */
    KAIJU_KAIJU2KRONA( ch_input_for_kaiju2krona.profiles, ch_input_for_kaiju2krona.db )
    ch_krona_text = ch_krona_text.mix( KAIJU_KAIJU2KRONA.out.txt )
    ch_versions = ch_versions.mix( KAIJU_KAIJU2KRONA.out.versions.first() )

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
        .map{[[id: it[0]['db_name'], tool: it[0]['tool']], it[1]]}
        .groupTuple()

    KRONA_KTIMPORTTEXT( ch_krona_text_for_import )
    ch_krona_html = ch_krona_html.mix( KRONA_KTIMPORTTEXT.out.html )
    ch_versions = ch_versions.mix( KRONA_KTIMPORTTEXT.out.versions.first() )

    /*
        Convert MALT/MEGAN RMA2INFO files into html Krona visualisations
    */
    if ( params.krona_taxonomy_directory ) {
        MEGAN_RMA2INFO_KRONA ( ch_input_classifications.malt, false )
        GUNZIP ( MEGAN_RMA2INFO_KRONA.out.txt )
        ch_krona_taxonomy_for_input = GUNZIP.out.gunzip
            .map{[[id: it[0]['db_name'], tool: it[0]['tool']], it[1]]}
            .groupTuple()

        KRONA_KTIMPORTTAXONOMY ( ch_krona_taxonomy_for_input, file(params.krona_taxonomy_directory, checkExists: true) )
        ch_krona_html.mix( KRONA_KTIMPORTTAXONOMY.out.html )
        ch_versions = ch_versions.mix( MEGAN_RMA2INFO_KRONA.out.versions.first() )
        ch_versions = ch_versions.mix( KRONA_KTIMPORTTAXONOMY.out.versions.first() )
    }

    emit:
    html = ch_krona_html
    versions = ch_versions
}
