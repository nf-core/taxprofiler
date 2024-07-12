//
// Process long raw reads with porechop
//

include { FASTQC as FASTQC_PROCESSED       } from '../../modules/nf-core/fastqc/main'
include { FALCO as FALCO_PROCESSED         } from '../../modules/nf-core/falco/main'

include { LONGREAD_ADAPTERREMOVAL          } from './longread_adapterremoval.nf'
include { LONGREAD_FILTERING               } from './longread_filtering.nf'

workflow LONGREAD_PREPROCESSING {
    take:
    reads

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( !params.longread_qc_skipadaptertrim && params.longread_qc_skipqualityfilter) {
        ch_processed_reads = LONGREAD_ADAPTERREMOVAL ( reads ).reads
        ch_versions        = ch_versions.mix(LONGREAD_ADAPTERREMOVAL.out.versions.first())
        ch_multiqc_files   =  ch_multiqc_files.mix( LONGREAD_ADAPTERREMOVAL.out.mqc )
    } else if ( params.longread_qc_skipadaptertrim && !params.longread_qc_skipqualityfilter) {
        ch_processed_reads = LONGREAD_FILTERING ( reads ).reads
        ch_versions        = ch_versions.mix(LONGREAD_FILTERING.out.versions.first())
        ch_multiqc_files   =  ch_multiqc_files.mix( LONGREAD_FILTERING.out.mqc )
    } else {
        LONGREAD_ADAPTERREMOVAL ( reads )
        ch_clipped_reads   = LONGREAD_ADAPTERREMOVAL.out.reads
            .map { meta, reads -> [ meta + [single_end: true], reads ] }
        ch_processed_reads = LONGREAD_FILTERING ( ch_clipped_reads ).reads
        ch_versions        = ch_versions.mix(LONGREAD_ADAPTERREMOVAL.out.versions.first())
        ch_versions        = ch_versions.mix(LONGREAD_FILTERING.out.versions.first())
        ch_multiqc_files   =  ch_multiqc_files.mix( LONGREAD_ADAPTERREMOVAL.out.mqc )
        ch_multiqc_files   =  ch_multiqc_files.mix( LONGREAD_FILTERING.out.mqc )
    }

    if (params.preprocessing_qc_tool == 'fastqc') {
        FASTQC_PROCESSED ( ch_processed_reads )
        ch_versions      = ch_versions.mix( FASTQC_PROCESSED.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip )

    } else if (params.preprocessing_qc_tool == 'falco') {
        FALCO_PROCESSED ( ch_processed_reads )
        ch_versions      = ch_versions.mix( FALCO_PROCESSED.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FALCO_PROCESSED.out.txt )
    }

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
