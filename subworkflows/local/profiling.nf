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
include { METAPHLAN3_METAPHLAN3                         } from '../../modules/nf-core/metaphlan3/metaphlan3/main'
include { KAIJU_KAIJU                                   } from '../../modules/nf-core/kaiju/kaiju/main'
include { KAIJU_KAIJU2TABLE as KAIJU_KAIJU2TABLE_SINGLE } from '../../modules/nf-core/kaiju/kaiju2table/main'
include { DIAMOND_BLASTX                                } from '../../modules/nf-core/diamond/blastx/main'
include { MOTUS_PROFILE                                 } from '../../modules/nf-core/motus/profile/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ                } from '../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'

workflow PROFILING {
    take:
    reads // [ [ meta ], [ reads ] ]
    databases // [ [ meta ], path ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_raw_classifications  = Channel.empty()
    ch_raw_profiles         = Channel.empty()

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
                metaphlan3: it[2]['tool'] == 'metaphlan3'
                motus: it[2]['tool'] == 'motus'
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
                    def temp_meta = [ id: meta['db_name'] ]

                    // Extend database parameters to specify whether to save alignments or not
                    def new_db_meta = db_meta.clone()
                    def sam_format = params.malt_save_reads ? ' --alignments ./ -za false' : ""
                    new_db_meta['db_params'] = db_meta['db_params'] + sam_format

                    // Combine reduced sample metadata with updated database parameters metadata,
                    // make sure id is db_name for publishing purposes.
                    def new_meta = temp_meta + new_db_meta
                    new_meta['id'] = new_meta['db_name']

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
                                            def meta_new = meta.clone()
                                            meta_new['db_name'] = meta.id
                                            meta_new['id'] = rma.baseName
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
                                        def db_meta_new = db_meta.clone()

                                        // Only take second element if one exists
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
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

        KRAKEN2_KRAKEN2 ( ch_input_for_kraken2.reads, ch_input_for_kraken2.db, params.kraken2_save_reads, params.kraken2_save_readclassification )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKEN2_KRAKEN2.out.classified_reads_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix(
            KRAKEN2_KRAKEN2.out.report
                // Set the tool to be strictly 'kraken2' instead of potentially 'bracken' for downstream use.
                // Will remain distinct from 'pure' Kraken2 results due to distinct database names in file names.
                .map { meta, report -> [meta + [tool: 'kraken2'], report]}
        )

    }

    if ( params.run_kraken2 && params.run_bracken ) {
        // Remove files from 'pure' kraken2 runs, so only those aligned against Bracken & kraken2 database are used.
        def ch_kraken2_output = KRAKEN2_KRAKEN2.out.report
            .filter {
                meta, report ->
                    if ( meta['instrument_platform'] == 'OXFORD_NANOPORE' ) log.warn "[nf-core/taxprofiler] Bracken has not been evaluated for Nanopore data. Skipping Bracken for sample ${meta.id}."
                    meta['tool'] == 'bracken' && meta['instrument_platform'] != 'OXFORD_NANOPORE'
            }

        // If necessary, convert the eight column output to six column output.
        if (params.kraken2_save_minimizers) {
            ch_kraken2_output = KRAKEN2_STANDARD_REPORT(ch_kraken2_output).report
        }

        // Extract the database name to combine by.
        ch_bracken_databases = databases
            .filter { meta, db -> meta['tool'] == 'bracken' }
            .map { meta, db -> [meta['db_name'], meta, db] }

        // Combine back with the reads
        ch_input_for_bracken = ch_kraken2_output
            .map { meta, report -> [meta['db_name'], meta, report] }
            .combine(ch_bracken_databases, by: 0)
            .map {

                key, meta, reads, db_meta, db ->
                    def db_meta_new = db_meta.clone()

                    // Have to pick second element if using bracken, as first element
                    // contains kraken parameters
                    if ( db_meta['tool'] == 'bracken' ) {

                        // Only take second element if one exists
                        def parsed_params = db_meta_new['db_params'].split(";")
                        if ( parsed_params.size() == 2 ) {
                            db_meta_new['db_params'] =  parsed_params[1]
                        } else {
                            db_meta_new['db_params'] = ""
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

        CENTRIFUGE_CENTRIFUGE ( ch_input_for_centrifuge.reads, ch_input_for_centrifuge.db, params.centrifuge_save_reads, params.centrifuge_save_reads, params.centrifuge_save_reads  )
        CENTRIFUGE_KREPORT (CENTRIFUGE_CENTRIFUGE.out.report, ch_input_for_centrifuge.db)
        ch_versions            = ch_versions.mix( CENTRIFUGE_CENTRIFUGE.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( CENTRIFUGE_CENTRIFUGE.out.results )
        ch_raw_profiles        = ch_raw_profiles.mix( CENTRIFUGE_KREPORT.out.kreport )
        ch_multiqc_files       = ch_multiqc_files.mix( CENTRIFUGE_KREPORT.out.kreport )

    }

    if ( params.run_metaphlan3 ) {

        ch_input_for_metaphlan3 = ch_input_for_profiling.metaphlan3
                            .filter{
                                if (it[0].is_fasta) log.warn "[nf-core/taxprofiler] MetaPhlAn3 currently does not accept FASTA files as input. Skipping MetaPhlAn3 for sample ${it[0].id}."
                                !it[0].is_fasta
                            }
                            .multiMap {
                                it ->
                                    reads: [it[0] + it[2], it[1]]
                                    db: it[3]
                            }

        METAPHLAN3_METAPHLAN3 ( ch_input_for_metaphlan3.reads, ch_input_for_metaphlan3.db )
        ch_versions        = ch_versions.mix( METAPHLAN3_METAPHLAN3.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( METAPHLAN3_METAPHLAN3.out.profile )

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

        KAIJU_KAIJU2TABLE_SINGLE ( KAIJU_KAIJU.out.results, ch_input_for_kaiju.db, params.kaiju_taxon_rank)
        ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE_SINGLE.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE_SINGLE.out.summary )
        ch_raw_profiles    = ch_raw_profiles.mix( KAIJU_KAIJU2TABLE_SINGLE.out.summary )
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
                                    .multiMap {
                                        single_meta, reads, db_meta, db ->
                                            reads: [ single_meta + db_meta, reads.flatten() ]
                                            db: db
                                }
        // Hardcode to _always_ produce the report file (which is our basic output, and goes into)
        KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( ch_input_for_krakenuniq.reads, ch_input_for_krakenuniq.db, params.krakenuniq_ram_chunk_size, params.krakenuniq_save_reads, true, params.krakenuniq_save_readclassifications )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
        ch_versions            = ch_versions.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.classified_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )

    }

    emit:
    classifications = ch_raw_classifications
    profiles        = ch_raw_profiles    // channel: [ val(meta), [ reads ] ] - should be text files or biom
    versions        = ch_versions          // channel: [ versions.yml ]
    motus_version   = params.run_motus ? MOTUS_PROFILE.out.versions.first() : Channel.empty()
    mqc             = ch_multiqc_files
}
