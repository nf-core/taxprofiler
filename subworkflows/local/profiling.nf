//
// Run profiling
//

include { MALT_RUN                    } from '../../modules/nf-core/modules/malt/run/main'
include { KRAKEN2_KRAKEN2             } from '../../modules/nf-core/modules/kraken2/kraken2/main'
include { CENTRIFUGE_CENTRIFUGE       } from '../../modules/nf-core/modules/centrifuge/centrifuge/main'
include { METAPHLAN3                  } from '../../modules/nf-core/modules/metaphlan3/main'

workflow PROFILING {
    take:
    reads // [ [ meta ], [ reads ] ]
    databases // [ [ meta ], path ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

/*
        COMBINE READS WITH POSSIBLE DATABASES
    */

    // e.g. output [DUMP: reads_plus_db] [['id':'2612', 'run_accession':'combined', 'instrument_platform':'ILLUMINA', 'single_end':1], <reads_path>/2612.merged.fastq.gz, ['tool':'malt', 'db_name':'mal95', 'db_params':'"-id 90"'], <db_path>/malt90]
    ch_input_for_profiling = reads
            .combine(databases)
            .map {
                meta, reads ->
                    def meta_new = meta.clone()
                        pairtype = meta_new['single_end'] ? '_se' : '_pe'
                        meta_new['id'] =  meta_new['id'] + pairtype
                    [meta_new, reads]
            }
            .branch {
                malt:    it[2]['tool'] == 'malt'
                kraken2: it[2]['tool'] == 'kraken2'
                metaphlan3: it[2]['tool'] == 'metaphlan3'
                centrifuge: it[2]['tool'] == 'centrifuge'
                unknown: true
            }

    /*
        PREPARE PROFILER INPUT CHANNELS
    */

    // Each tool as a slightly different input structure and generally separate
    // input channels for reads vs databases. We restructure the channel tuple
    // for each tool and make liberal use of multiMap to keep reads/databases
    // channel element order in sync with each other

    // MALT: We groupTuple to have all samples in one channel for MALT as database
    // loading takes a long time, so we only want to run it once per database
    // TODO document somewhere we only accept illumina short reads for MALT?
    ch_input_for_malt =  ch_input_for_profiling.malt
                            .filter { it[0]['instrument_platform'] == 'ILLUMINA' }
                            .map {
                                it ->
                                    def temp_meta =  [ id: it[2]['db_name']]  + it[2]
                                    def db = it[3]
                                    [ temp_meta, it[1], db ]
                            }
                            .groupTuple(by: [0,2])
                            .multiMap {
                                it ->
                                    reads: [ it[0], it[1].flatten() ]
                                    db: it[2]
                            }

    // All subsequent tools can easily run on a per-sample basis

    ch_input_for_kraken2 =  ch_input_for_profiling.kraken2
                            .multiMap {
                                it ->
                                    reads: [ it[0] + it[2], it[1] ]
                                    db: it[3]
                            }

    ch_input_for_centrifuge =  ch_input_for_profiling.centrifuge
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

    ch_input_for_metaphlan3 = ch_input_for_profiling.metaphlan3
                            .multiMap {
                                it ->
                                    reads: [it[0] + it[2], it[1]]
                                    db: it[3]
                            }

    /*
        RUN PROFILING
    */

    if ( params.run_malt ) {
        MALT_RUN ( ch_input_for_malt.reads, params.malt_mode, ch_input_for_malt.db )
        ch_multiqc_files = ch_multiqc_files.mix( MALT_RUN.out.log.collect{it[1]}.ifEmpty([])  )
        ch_versions = ch_versions.mix( MALT_RUN.out.versions.first() )
    }

    if ( params.run_kraken2 ) {
        KRAKEN2_KRAKEN2 ( ch_input_for_kraken2.reads, ch_input_for_kraken2.db  )
        ch_multiqc_files = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.txt.collect{it[1]}.ifEmpty([])  )
        ch_versions = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
    }

    if ( params.run_centrifuge ) {
        CENTRIFUGE_CENTRIFUGE ( ch_input_for_centrifuge.reads, ch_input_for_centrifuge.db, params.centrifuge_save_unaligned, params.centrifuge_save_aligned, params.centrifuge_sam_format  )
        ch_versions = ch_versions.mix( CENTRIFUGE_CENTRIFUGE.out.versions.first() )
    }

    if ( params.run_metaphlan3 ) {
        METAPHLAN3 ( ch_input_for_metaphlan3.reads, ch_input_for_metaphlan3.db )
        ch_versions = ch_versions.mix( METAPHLAN3.out.versions.first() )
    }


    emit:
    // TODO work out if there is enough standardisation of output to export as one?
    //output    = ch_filtered_reads    // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

