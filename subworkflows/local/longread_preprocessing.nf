//
// Process long raw reads with porechop
//

include { FASTQC as FASTQC_PROCESSED } from '../../modules/nf-core/fastqc/main'
include { FALCO as FALCO_PROCESSED   } from '../../modules/nf-core/falco/main'

include { PORECHOP_PORECHOP          } from '../../modules/nf-core/porechop/porechop/main'
include { FILTLONG                   } from '../../modules/nf-core/filtlong/main'

workflow LONGREAD_PREPROCESSING {
    take:
    reads

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( !params.longread_qc_skipadaptertrim && params.longread_qc_skipqualityfilter) {
        PORECHOP_PORECHOP ( reads )

        ch_processed_reads = PORECHOP_PORECHOP.out.reads
            .map { meta, reads -> [ meta + [single_end: 1], reads ] }

        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_PORECHOP.out.log )

    } else if ( params.longread_qc_skipadaptertrim && !params.longread_qc_skipqualityfilter) {

        ch_processed_reads = FILTLONG ( reads.map { meta, reads -> [meta, [], reads ] } )
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( FILTLONG.out.log )

    } else {
        PORECHOP_PORECHOP ( reads )
        ch_clipped_reads = PORECHOP_PORECHOP.out.reads
            .map { meta, reads -> [ meta + [single_end: 1], reads ] }

        ch_processed_reads = FILTLONG ( ch_clipped_reads.map { meta, reads -> [ meta, [], reads ] } ).reads

        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_PORECHOP.out.log )
        ch_multiqc_files = ch_multiqc_files.mix( FILTLONG.out.log )
    }

    if (params.preprocessing_qc_tool == 'fastqc') {
        FASTQC_PROCESSED ( ch_processed_reads )
        ch_versions = ch_versions.mix( FASTQC_PROCESSED.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip )

    } else if (params.preprocessing_qc_tool == 'falco') {
        FALCO_PROCESSED ( ch_processed_reads )
        ch_versions = ch_versions.mix( FALCO_PROCESSED.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FALCO_PROCESSED.out.txt )
    }

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

