//
// Run profiling
//

include { MALT_RUN                              } from '../../modules/nf-core/modules/malt/run/main'
include { MEGAN_RMA2INFO as MEGAN_RMA2INFO_TSV  } from '../../modules/nf-core/modules/megan/rma2info/main'
include { KRAKEN2_KRAKEN2                       } from '../../modules/nf-core/modules/kraken2/kraken2/main'
include { CENTRIFUGE_CENTRIFUGE                 } from '../../modules/nf-core/modules/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT                    } from '../../modules/nf-core/modules/centrifuge/kreport/main'
include { METAPHLAN3                            } from '../../modules/nf-core/modules/metaphlan3/metaphlan3/main'
include { KAIJU_KAIJU                           } from '../../modules/nf-core/modules/kaiju/kaiju/main'
include { DIAMOND_BLASTX                        } from '../../modules/nf-core/modules/diamond/blastx/main'
include { MOTUS_PROFILE                         } from '../../modules/nf-core/modules/motus/profile/main'

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
                    def meta_new = meta.clone()
                        pairtype = meta_new['single_end'] ? '_se' : '_pe'
                        meta_new['id'] =  meta_new['id'] + pairtype
                    [meta_new, reads]
            }
            .combine(databases)
            .branch {
                malt:    it[2]['tool'] == 'malt'
                kraken2: it[2]['tool'] == 'kraken2'
                metaphlan3: it[2]['tool'] == 'metaphlan3'
                centrifuge: it[2]['tool'] == 'centrifuge'
                kaiju: it[2]['tool'] == 'kaiju'
                diamond: it[2]['tool'] == 'diamond'
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
                                .filter { it[0]['instrument_platform'] == 'ILLUMINA' }
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
                                    it ->
                                        reads: [ it[0], it[1].flatten() ]
                                        db: it[2]
                                }

        MALT_RUN ( ch_input_for_malt.reads, params.malt_mode, ch_input_for_malt.db )

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

    if ( params.run_kraken2 ) {

        ch_input_for_kraken2 =  ch_input_for_profiling.kraken2
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

        KRAKEN2_KRAKEN2 ( ch_input_for_kraken2.reads, ch_input_for_kraken2.db, params.kraken2_save_reads, params.kraken2_save_readclassification )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKEN2_KRAKEN2.out.classified_reads_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix( KRAKEN2_KRAKEN2.out.report )

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

        METAPHLAN3 ( ch_input_for_metaphlan3.reads, ch_input_for_metaphlan3.db )
        ch_versions        = ch_versions.mix( METAPHLAN3.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( METAPHLAN3.out.biom )

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

    emit:
    classifications = ch_raw_classifications
    profiles        = ch_raw_profiles    // channel: [ val(meta), [ reads ] ] - should be text files or biom
    versions        = ch_versions          // channel: [ versions.yml ]
    motus_version   = params.run_motus ? MOTUS_PROFILE.out.versions.first() : Channel.empty()
    mqc             = ch_multiqc_files
}
