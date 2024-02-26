//
// Run profiling
//

include { MALT_RUN                                      } from '../../modules/nf-core/malt/run/main'
include { MEGAN_RMA2INFO as MEGAN_RMA2INFO_TSV          } from '../../modules/nf-core/megan/rma2info/main'
include { KRAKEN2_KRAKEN2                               } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_STANDARD_REPORT                       } from '../../modules/local/kraken2_standard_report'
include { BRACKEN_BRACKEN                               } from '../../modules/nf-core/bracken/bracken/main'
include { CENTRIFUGE_CENTRIFUGE                         } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT                            } from '../../modules/nf-core/centrifuge/kreport/main'
include { METAPHLAN_METAPHLAN                           } from '../../modules/nf-core/metaphlan/metaphlan/main'
include { KAIJU_KAIJU                                   } from '../../modules/nf-core/kaiju/kaiju/main'
include { KAIJU_KAIJU2TABLE as KAIJU_KAIJU2TABLE_SINGLE } from '../../modules/nf-core/kaiju/kaiju2table/main'
include { DIAMOND_BLASTX                                } from '../../modules/nf-core/diamond/blastx/main'
include { MOTUS_PROFILE                                 } from '../../modules/nf-core/motus/profile/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ                } from '../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'
include { KMCP_SEARCH                                   } from '../../modules/nf-core/kmcp/search/main'
include { KMCP_PROFILE                                  } from '../../modules/nf-core/kmcp/profile/main'
include { GANON_CLASSIFY                                } from '../../modules/nf-core/ganon/classify/main'
include { GANON_REPORT                                  } from '../../modules/nf-core/ganon/report/main'


// Custom Functions

/**
* Combine profiles with their original database, then separate into two channels.
*
* The channel elements are assumed to be tuples one of [ meta, profile ], and the
* database to be of [db_key, meta, database_file].
*
* @param ch_profile A channel containing a meta and the profilign report of a given profiler
* @param ch_database A channel containing a key, the database meta, and the database file/folders itself
* @return A multiMap'ed output channel with two sub channels, one with the profile and the other with the db
*/
def combineProfilesWithDatabase(ch_profile, ch_database) {

return ch_profile
    .map { meta, profile -> [meta.db_name, meta, profile] }
    .combine(ch_database, by: 0)
    .multiMap {
        key, meta, profile, db_meta, db ->
            profile: [meta, profile]
            db: db
    }
}

workflow PROFILING {
    take:
    reads // [ [ meta ], [ reads ] ]
    databases // [ [ meta ], path ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_raw_classifications  = Channel.empty() // These per-read ID taxonomic assingment
    ch_raw_profiles         = Channel.empty() // These are count tables

/*
        COMBINE READS WITH POSSIBLE DATABASES
    */

    // e.g. output [DUMP: reads_plus_db] [['id':'2612', 'run_accession':'combined', 'instrument_platform':'ILLUMINA', 'single_end':1], <reads_path>/2612.merged.fastq.gz, ['tool':'malt', 'db_name':'mal95', 'db_params':'"-id 90"'], <db_path>/malt90]
    ch_input_for_profiling = reads
            .map {
                meta, reads ->
                    [meta + [id: "${meta.id}${meta.single_end ? '_se' : '_pe'}"], reads]
            }
            .combine(databases)
            .branch {
                centrifuge: it[2]['tool'] == 'centrifuge'
                diamond: it[2]['tool'] == 'diamond'
                kaiju: it[2]['tool'] == 'kaiju'
                kraken2: it[2]['tool'] == 'kraken2' || it[2]['tool'] == 'bracken' // to reuse the kraken module to produce the input data for bracken
                krakenuniq: it[2]['tool'] == 'krakenuniq'
                malt:    it[2]['tool'] == 'malt'
                metaphlan: it[2]['tool'] == 'metaphlan'
                motus: it[2]['tool'] == 'motus'
                kmcp: it[2]['tool'] == 'kmcp'
                ganon: it[2]['tool'] == 'ganon'
                unknown: true
            }

    /*
        PREPARE PROFILER INPUT CHANNELS & RUN PROFILING
    */

    // Each tool as a slightly different input structure and generally separate
    // input channels for reads vs databases. We restructure the channel tuple
    // for each tool and make liberal use of multiMap to keep reads/databases
    // channel element order in sync with each other

    if ( params.run_malt ) {

        // MALT: We groupTuple to have all samples in one channel for MALT as database
        // loading takes a long time, so we only want to run it once per database
        ch_input_for_malt =  ch_input_for_profiling.malt
            .map {
                meta, reads, db_meta, db ->

                    // Reset entire input meta for MALT to just database name,
                    // as we don't run run on a per-sample basis due to huge datbaases
                    // so all samples are in one run and so sample-specific metadata
                    // unnecessary. Set as database name to prevent `null` job ID and prefix.
                    def new_meta = db_meta + [ id: db_meta.db_name ]

                    // Extend database parameters to specify whether to save alignments or not
                    def sam_format = params.malt_save_reads ? ' --alignments ./ -za false' : ""
                    new_meta.db_params = db_meta.db_params + sam_format

                    [ new_meta, reads, db ]

            }
            .groupTuple(by: [0,2])
            .multiMap {
                meta, reads, db ->
                    reads: [ meta, reads.flatten() ]
                    db: db
            }

        MALT_RUN ( ch_input_for_malt.reads, ch_input_for_malt.db )

        ch_maltrun_for_megan = MALT_RUN.out.rma6
                                .transpose()
                                .map{
                                    meta, rma ->
                                            // re-extract meta from file names, use filename without rma to
                                            // ensure we keep paired-end information in downstream filenames
                                            // when no pair-merging
                                            def meta_new = meta + [db_name: meta.id, id: rma.baseName]

                                        [ meta_new, rma ]
                                }

        MEGAN_RMA2INFO_TSV (ch_maltrun_for_megan, params.malt_generate_megansummary )
        ch_multiqc_files       = ch_multiqc_files.mix( MALT_RUN.out.log )
        ch_versions            = ch_versions.mix( MALT_RUN.out.versions.first(), MEGAN_RMA2INFO_TSV.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( ch_maltrun_for_megan )
        ch_raw_profiles        = ch_raw_profiles.mix( MEGAN_RMA2INFO_TSV.out.txt )

    }

    if ( params.run_kraken2 || params.run_bracken ) {
        // Have to pick first element of db_params if using bracken,
        // as db sheet for bracken must have ; sep list to
        // distinguish between kraken and bracken parameters
        ch_input_for_kraken2 = ch_input_for_profiling.kraken2
                                .map {
                                    meta, reads, db_meta, db ->

                                        // Only take first element if one exists
                                        def parsed_params = db_meta['db_params'].split(";")
                                        if ( parsed_params.size() == 2 ) {
                                            db_meta_new = db_meta + [db_params: parsed_params[0]]
                                        } else if ( parsed_params.size() == 0 ) {
                                            db_meta_new = db_meta + [db_params: ""]
                                        } else {
                                            db_meta_new = db_meta + [db_params: parsed_params[0]]
                                        }

                                    [ meta, reads, db_meta_new, db ]
                                }
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

        KRAKEN2_KRAKEN2 ( ch_input_for_kraken2.reads, ch_input_for_kraken2.db, params.kraken2_save_reads, params.kraken2_save_readclassifications )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKEN2_KRAKEN2.out.classified_reads_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix(
            KRAKEN2_KRAKEN2.out.report
                // Rename tool in the meta for the for-bracken files to disambiguate from only-kraken2 results in downstream steps.
                // Note may need to rename back to to just bracken in those downstream steps depending on context.
                .map {
                    meta, report ->
                        def new_tool =
                    [meta + [tool: meta.tool == 'bracken' ? 'kraken2-bracken' : meta.tool], report]
                }
        )

    }

    if ( params.run_kraken2 && params.run_bracken ) {
        // Remove files from 'pure' kraken2 runs, so only those aligned against Bracken & kraken2 database are used.
        def ch_kraken2_output = KRAKEN2_KRAKEN2.out.report
            .filter {
                meta, report ->
                    if ( meta.instrument_platform == 'OXFORD_NANOPORE' ) log.warn "[nf-core/taxprofiler] Bracken has not been evaluated for Nanopore data. Skipping Bracken for sample ${meta.id}."
                    meta.tool == 'bracken' && meta.instrument_platform != 'OXFORD_NANOPORE'
            }

        // If necessary, convert the eight column output to six column output.
        if (params.kraken2_save_minimizers) {
            ch_kraken2_output = KRAKEN2_STANDARD_REPORT(ch_kraken2_output).report
        }

        // Extract the database name to combine by.
        ch_bracken_databases = databases
            .filter { meta, db -> meta.tool == 'bracken' }
            .map { meta, db -> [meta.db_name, meta, db] }

        // Combine back with the reads
        ch_input_for_bracken = ch_kraken2_output
            .map { meta, report -> [meta.db_name, meta, report] }
            .combine(ch_bracken_databases, by: 0)
            .map {

                key, meta, reads, db_meta, db ->

                    // // Have to make a completely fresh copy here as otherwise
                    // // was getting db_param loss due to upstream meta parsing at
                    // // kraken2 input channel manipulation step
                    def db_meta_keys = db_meta.keySet()
                    def db_meta_new = db_meta.subMap(db_meta_keys)

                    // Have to pick second element if using bracken, as first element
                    // contains kraken parameters
                    if ( db_meta.tool == 'bracken' ) {

                        // Only take second element if one exists
                        def parsed_params = db_meta['db_params'].split(";")

                        if ( parsed_params.size() == 2 ) {
                            db_meta_new = db_meta + [ db_params: parsed_params[1] ]
                        } else {
                            db_meta_new = db_meta + [ db_params: "" ]
                        }

                    } else {
                        db_meta_new['db_params']
                    }

                [ key, meta, reads, db_meta_new, db ]
            }
            .multiMap { key, meta, report, db_meta, db ->
                report: [meta + db_meta, report]
                db: db
            }

        BRACKEN_BRACKEN(ch_input_for_bracken.report, ch_input_for_bracken.db)
        ch_versions     = ch_versions.mix(BRACKEN_BRACKEN.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(BRACKEN_BRACKEN.out.reports)

    }

    if ( params.run_centrifuge ) {

        ch_input_for_centrifuge =  ch_input_for_profiling.centrifuge
                                .filter{
                                    if (it[0].is_fasta) log.warn "[nf-core/taxprofiler] Centrifuge currently does not accept FASTA files as input. Skipping Centrifuge for sample ${it[0].id}."
                                    !it[0].is_fasta
                                }
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

        CENTRIFUGE_CENTRIFUGE ( ch_input_for_centrifuge.reads, ch_input_for_centrifuge.db, params.centrifuge_save_reads, params.centrifuge_save_reads  )
        ch_versions            = ch_versions.mix( CENTRIFUGE_CENTRIFUGE.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( CENTRIFUGE_CENTRIFUGE.out.results )

        // Ensure the correct database goes with the generated report for KREPORT
        ch_database_for_centrifugekreport = databases
                                                .filter { meta, db -> meta.tool == 'centrifuge' }
                                                .map { meta, db -> [meta.db_name, meta, db] }

        ch_input_for_centrifuge_kreport = combineProfilesWithDatabase(CENTRIFUGE_CENTRIFUGE.out.results, ch_database_for_centrifugekreport)

        // Generate profile
        CENTRIFUGE_KREPORT (ch_input_for_centrifuge_kreport.profile, ch_input_for_centrifuge_kreport.db)
        ch_versions            = ch_versions.mix( CENTRIFUGE_KREPORT.out.versions.first() )
        ch_raw_profiles        = ch_raw_profiles.mix( CENTRIFUGE_KREPORT.out.kreport )
        ch_multiqc_files       = ch_multiqc_files.mix( CENTRIFUGE_KREPORT.out.kreport )

    }

    if ( params.run_metaphlan ) {

        ch_input_for_metaphlan = ch_input_for_profiling.metaphlan
                            .multiMap {
                                it ->
                                    reads: [it[0] + it[2], it[1]]
                                    db: it[3]
                            }

        METAPHLAN_METAPHLAN ( ch_input_for_metaphlan.reads, ch_input_for_metaphlan.db )
        ch_versions        = ch_versions.mix( METAPHLAN_METAPHLAN.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( METAPHLAN_METAPHLAN.out.profile )

    }

    if ( params.run_kaiju ) {

        ch_input_for_kaiju = ch_input_for_profiling.kaiju
                            .multiMap {
                                it ->
                                    reads: [it[0] + it[2], it[1]]
                                    db: it[3]
                            }

        KAIJU_KAIJU ( ch_input_for_kaiju.reads, ch_input_for_kaiju.db)
        ch_versions = ch_versions.mix( KAIJU_KAIJU.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KAIJU_KAIJU.out.results )

        // Ensure the correct database goes with the generated report for KAIJU2TABLE
        ch_database_for_kaiju2table = databases
                                                .filter { meta, db -> meta.tool == 'kaiju' }
                                                .map { meta, db -> [meta.db_name, meta, db] }

        ch_input_for_kaiju2table = combineProfilesWithDatabase(KAIJU_KAIJU.out.results, ch_database_for_kaiju2table)
        // Generate profile
        KAIJU_KAIJU2TABLE_SINGLE ( ch_input_for_kaiju2table.profile, ch_input_for_kaiju2table.db, params.kaiju_taxon_rank)
        ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE_SINGLE.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE_SINGLE.out.summary )
        ch_raw_profiles  = ch_raw_profiles.mix( KAIJU_KAIJU2TABLE_SINGLE.out.summary )
    }

    if ( params.run_diamond ) {

        ch_input_for_diamond = ch_input_for_profiling.diamond
                                .multiMap {
                                    it ->
                                        reads: [it[0] + it[2], it[1]]
                                        db: it[3]
                                }

        // diamond only accepts single output file specification, therefore
        // this will replace output file!
        ch_diamond_reads_format = params.diamond_save_reads ? 'sam' : params.diamond_output_format

        DIAMOND_BLASTX ( ch_input_for_diamond.reads, ch_input_for_diamond.db, ch_diamond_reads_format , [] )
        ch_versions        = ch_versions.mix( DIAMOND_BLASTX.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( DIAMOND_BLASTX.out.tsv )
        ch_multiqc_files   = ch_multiqc_files.mix( DIAMOND_BLASTX.out.log )

    }

    if ( params.run_motus ) {

        ch_input_for_motus = ch_input_for_profiling.motus
                                .filter{
                                    if (it[0].is_fasta) log.warn "[nf-core/taxprofiler] mOTUs currently does not accept FASTA files as input. Skipping mOTUs for sample ${it[0].id}."
                                    !it[0].is_fasta
                                }
                                .multiMap {
                                    it ->
                                        reads: [it[0] + it[2], it[1]]
                                        db: it[3]
                                }

        MOTUS_PROFILE ( ch_input_for_motus.reads, ch_input_for_motus.db )
        ch_versions        = ch_versions.mix( MOTUS_PROFILE.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( MOTUS_PROFILE.out.out )
        ch_multiqc_files   = ch_multiqc_files.mix( MOTUS_PROFILE.out.log )
    }

    if ( params.run_krakenuniq ) {
        ch_input_for_krakenuniq =  ch_input_for_profiling.krakenuniq
            .map {
                meta, reads, db_meta, db ->
                    [[id: db_meta.db_name, single_end: meta.single_end], reads, db_meta, db]
            }
            .groupTuple(by: [0,2,3])
            .flatMap { single_meta, reads, db_meta, db ->
                def batches = reads.collate(params.krakenuniq_batch_size)
                return batches.collect { batch -> [ single_meta + db_meta, batch.flatten(), db ]}
            }
            .multiMap {
                meta, reads, db ->
                    reads: [ meta, reads ]
                    db: db
            }
        // Hardcode to _always_ produce the report file (which is our basic output, and goes into)
        KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( ch_input_for_krakenuniq.reads, ch_input_for_krakenuniq.db, params.krakenuniq_ram_chunk_size, params.krakenuniq_save_reads, true, params.krakenuniq_save_readclassifications )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
        ch_versions            = ch_versions.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.classified_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )

    }

    if (params.run_kmcp) {

            ch_input_for_kmcp = ch_input_for_profiling.kmcp
                                .filter {
                                    meta, reads, meta_db, db ->
                                        if ( meta['instrument_platform'] == 'OXFORD_NANOPORE' ) log.warn "[nf-core/taxprofiler] KMCP is only suitable for short-read metagenomic profiling, with much lower sensitivity on long-read datasets. Skipping KMCP for sample ${meta.id}."
                                        meta_db['tool'] == 'kmcp' && meta['instrument_platform'] != 'OXFORD_NANOPORE'
                                    }
                                .map {
                                    meta, reads, db_meta, db ->
                                        def db_meta_keys = db_meta.keySet()
                                        def db_meta_new = db_meta.subMap(db_meta_keys)

                                        // Split the string, the arguments before semicolon should be parsed into kmcp search
                                        def parsed_params = db_meta_new['db_params'].split(";")
                                        if ( parsed_params.size() == 2 ) {
                                            db_meta_new['db_params'] = parsed_params[0]
                                        } else if ( parsed_params.size() == 0 ) {
                                            db_meta_new['db_params'] = ""
                                        } else {
                                            db_meta_new['db_params'] = parsed_params[0]
                                        }

                                    [ meta, reads, db_meta_new, db ]
                                }
                                .multiMap  {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

            KMCP_SEARCH ( ch_input_for_kmcp.db, ch_input_for_kmcp.reads )

            ch_versions            = ch_versions.mix( KMCP_SEARCH.out.versions.first() )
            ch_raw_classifications = ch_raw_classifications.mix(KMCP_SEARCH.out.result)

            ch_database_for_kmcp_profile = databases
                                                .filter { meta, db -> meta.tool == 'kmcp' }
                                                .map { meta, db -> [meta.db_name, meta, db] }

            ch_input_for_kmcp_profile = KMCP_SEARCH.out.result
                .map { meta, report -> [meta.db_name, meta, report] }
                .combine(ch_database_for_kmcp_profile, by: 0)
                .map {

                    key, meta, reads, db_meta, db ->

                        // Same as kraken2/bracken logic here. Arguments after semicolon are going into KMCP_PROFILE
                        def db_meta_keys = db_meta.keySet()
                        def db_meta_new = db_meta.subMap(db_meta_keys)

                        def parsed_params = db_meta['db_params'].split(";")

                            if ( parsed_params.size() == 2 ) {
                                db_meta_new = db_meta + [ db_params: parsed_params[1] ]
                            } else {
                                db_meta_new = db_meta + [ db_params: "" ]
                            }

                    [ key, meta, reads, db_meta_new, db ]

            }
            .multiMap { key, meta, report, db_meta, db ->
                report: [meta + db_meta, report]
                db: db
            }

            //Generate kmcp profile
            KMCP_PROFILE( ch_input_for_kmcp_profile.report, ch_input_for_kmcp.db, params.kmcp_mode )
            ch_versions = ch_versions.mix( KMCP_PROFILE.out.versions.first() )
            ch_raw_profiles    = ch_raw_profiles.mix( KMCP_PROFILE.out.profile )
            ch_multiqc_files   = ch_multiqc_files.mix( KMCP_PROFILE.out.profile )
}


    if ( params.run_ganon ) {

        ch_input_for_ganonclassify = ch_input_for_profiling.ganon
                                .filter {
                                    meta, reads, meta_db, db ->
                                        if ( meta.instrument_platform == 'OXFORD_NANOPORE' ) log.warn "[nf-core/taxprofiler] Ganon has not been evaluated for Nanopore data. Skipping Ganon for sample ${meta.id}."
                                        meta_db.tool == 'ganon' && meta.instrument_platform != 'OXFORD_NANOPORE'
                                }
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

        ch_input_for_ganonclassify.reads

        GANON_CLASSIFY( ch_input_for_ganonclassify.reads, ch_input_for_ganonclassify.db )
        ch_versions = ch_versions.mix( GANON_CLASSIFY.out.versions.first() )

        ch_database_for_ganonreport = databases
                                        .filter { meta, db -> meta.tool == "ganon" }
                                        .map { meta, db -> [meta.db_name, meta, db] }

        ch_report_for_ganonreport = combineProfilesWithDatabase(GANON_CLASSIFY.out.report, ch_database_for_ganonreport)

        GANON_REPORT(ch_report_for_ganonreport.profile, ch_report_for_ganonreport.db)
        ch_versions            = ch_versions.mix( GANON_REPORT.out.versions.first() )

        // Might be flipped - check/define what is a profile vs raw classification
        ch_raw_profiles        = ch_raw_profiles.mix( GANON_REPORT.out.tre )
        ch_raw_classifications = ch_raw_classifications.mix( GANON_CLASSIFY.out.all )

    }

    emit:
    classifications = ch_raw_classifications
    profiles        = ch_raw_profiles    // channel: [ val(meta), [ reads ] ] - should be text files or biom
    versions        = ch_versions          // channel: [ versions.yml ]
    motus_version   = params.run_motus ? MOTUS_PROFILE.out.versions.first() : Channel.empty()
    mqc             = ch_multiqc_files
}
