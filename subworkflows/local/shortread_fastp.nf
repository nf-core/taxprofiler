/*
Process short raw reads with FastP
*/

include { FASTP as FASTP_SINGLE       } from '../../modules/nf-core/modules/fastp/main'
include { FASTP as FASTP_PAIRED       } from '../../modules/nf-core/modules/fastp/main'

workflow SHORTREAD_FASTP {
    take:
    reads // [[meta], [reads]]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files      = Channel.empty()

    ch_input_for_fastp = reads
                            .dump(tag: "pre-fastp_branch")
                            .branch{
                                single: it[0]['single_end'] == true
                                paired: it[0]['single_end'] == false
                            }

    ch_input_for_fastp.single.dump(tag: "input_fastp_single")
    ch_input_for_fastp.paired.dump(tag: "input_fastp_paired")

    FASTP_SINGLE ( ch_input_for_fastp.single, false, false )
    // Last parameter here turns on merging of PE data
    FASTP_PAIRED ( ch_input_for_fastp.paired, false, params.shortread_clipmerge_mergepairs )

    if ( params.shortread_clipmerge_mergepairs ) {
        // TODO update to replace meta suffix
        ch_fastp_reads_prepped = FASTP_PAIRED.out.reads_merged
                                    .mix( FASTP_SINGLE.out.reads )
                                    .map {
                                        meta, reads ->
                                            def meta_new = meta.clone()
                                            meta_new['single_end'] = 1
                                            [ meta_new, reads ]
                                        }
    } else {
        ch_fastp_reads_prepped = FASTP_PAIRED.out.reads
                                    .mix( FASTP_SINGLE.out.reads )
    }

    ch_versions = ch_versions.mix(FASTP_SINGLE.out.versions.first())
    ch_versions = ch_versions.mix(FASTP_PAIRED.out.versions.first())

    ch_processed_reads = ch_fastp_reads_prepped.dump(tag: "ch_fastp_reads_prepped")

    ch_multiqc_files = ch_multiqc_files.mix( FASTP_SINGLE.out.json.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( FASTP_PAIRED.out.json.collect{it[1]} )

    ch_multiqc_files.dump(tag: "preprocessing_fastp_mqc_final")

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

