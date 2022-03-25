/*
Process long raw reads with porechop
*/

include { FASTQC as FASTQC_PROCESSED } from '../../modules/nf-core/modules/fastqc/main'
include { PORECHOP                   } from '../../modules/nf-core/modules/porechop/main'

workflow LONGREAD_PREPROCESSING {
    take:
    reads

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    PORECHOP ( reads )

    ch_processed_reads = PORECHOP.out.reads
                                .dump(tag: "pre_fastqc_check")
                                .map {
                                        meta, reads ->
                                        def meta_new = meta.clone()
                                        meta_new['single_end'] = 1
                                        [ meta_new, reads ]
                                    }

    FASTQC_PROCESSED ( PORECHOP.out.reads )
    ch_versions = ch_versions.mix(PORECHOP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip.collect{it[1]} )


    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

