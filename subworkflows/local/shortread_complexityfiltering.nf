//
// Check input samplesheet and get read channels
//

include { BBMAP_BBDUK     } from '../../modules/nf-core/bbmap/bbduk/main'
include { PRINSEQPLUSPLUS } from '../../modules/nf-core/prinseqplusplus/main'

workflow SHORTREAD_COMPLEXITYFILTERING {
    take:
    reads // [ [ meta ], [ reads ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    // fastp complexity filtering is activated via modules.conf in shortread_preprocessing
    if ( params.shortread_complexityfilter_tool == 'bbduk' ) {
        ch_filtered_reads = BBMAP_BBDUK ( reads, [] ).reads
        ch_versions        =  ch_versions.mix( BBMAP_BBDUK.out.versions.first() )
        ch_multiqc_files   =  ch_multiqc_files.mix( BBMAP_BBDUK.out.log )
    } else if ( params.shortread_complexityfilter_tool == 'prinseqplusplus' ) {
        ch_filtered_reads = PRINSEQPLUSPLUS ( reads ).good_reads
        ch_versions        =  ch_versions.mix( PRINSEQPLUSPLUS.out.versions.first() )
    } else {
        ch_filtered_reads = reads
    }

    emit:
    reads    = ch_filtered_reads    // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

