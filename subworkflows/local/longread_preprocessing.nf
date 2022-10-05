//
// Process long raw reads with porechop
//

include { FASTQC as FASTQC_PROCESSED } from '../../modules/nf-core/modules/fastqc/main'
include { PORECHOP                   } from '../../modules/nf-core/modules/porechop/main'
include { FILTLONG                   } from '../../modules/nf-core/modules/filtlong/main'

workflow LONGREAD_PREPROCESSING {
    take:
    reads

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( !params.longread_qc_skipadaptertrim && params.longread_qc_skipqualityfilter) {
        PORECHOP ( reads )

        ch_processed_reads = PORECHOP.out.reads
                                        .map {
                                                meta, reads ->
                                                def meta_new = meta.clone()
                                                meta_new['single_end'] = 1
                                                [ meta_new, reads ]
                                        }

        ch_versions = ch_versions.mix(PORECHOP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP.out.log )

    } else if ( params.longread_qc_skipadaptertrim && !params.longread_qc_skipqualityfilter) {

        ch_processed_reads = FILTLONG ( reads.map{ meta, reads -> [meta, [], reads ]} )
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( FILTLONG.out.log )

    } else {
        PORECHOP ( reads )
        ch_clipped_reads = PORECHOP.out.reads
                                        .map {
                                                meta, reads ->
                                                def meta_new = meta.clone()
                                                meta_new['single_end'] = 1
                                                [ meta_new, reads ]
                                        }

        ch_processed_reads = FILTLONG ( ch_clipped_reads.map{ meta, reads -> [meta, [], reads ]} ).reads

        ch_versions = ch_versions.mix(PORECHOP.out.versions.first())
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP.out.log )
        ch_multiqc_files = ch_multiqc_files.mix( FILTLONG.out.log )
    }

    FASTQC_PROCESSED ( ch_processed_reads )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip )

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

