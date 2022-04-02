//
// Perform read trimming and merging
//


include { SHORTREAD_FASTP             } from './shortread_fastp'
include { SHORTREAD_ADAPTERREMOVAL    } from './shortread_adapterremoval'
include { FASTQC as FASTQC_PROCESSED       } from '../../modules/nf-core/modules/fastqc/main'

workflow SHORTREAD_PREPROCESSING {
    take:
    reads //  [ [ meta ], [ reads ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.shortread_clipmerge_tool == "fastp" ) {
        ch_processed_reads = SHORTREAD_FASTP ( reads ).reads
        ch_versions        =  ch_versions.mix( SHORTREAD_FASTP.out.versions )
        ch_multiqc_files   =  ch_multiqc_files.mix( SHORTREAD_FASTP.out.mqc )
    } else if ( params.shortread_clipmerge_tool == "adapterremoval" ) {
        ch_processed_reads = SHORTREAD_ADAPTERREMOVAL ( reads ).reads
        ch_versions        = ch_versions.mix( SHORTREAD_ADAPTERREMOVAL.out.versions )
        ch_multiqc_files   = ch_multiqc_files.mix( SHORTREAD_ADAPTERREMOVAL.out.mqc )
    } else {
        ch_processed_reads = reads
    }

    FASTQC_PROCESSED ( ch_processed_reads )
    ch_versions = ch_versions.mix( FASTQC_PROCESSED.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip )

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

