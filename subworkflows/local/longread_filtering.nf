//
// Perform filtering
//

include { FILTLONG                   } from '../../modules/nf-core/filtlong/main'
include { NANOQ                      } from '../../modules/nf-core/nanoq/main'

workflow LONGREAD_FILTERING {
    take:
    reads // [ [ meta ], [ reads ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    // fastp complexity filtering is activated via modules.conf in shortread_preprocessing
    if ( params.longread_filter_tool == 'filtlong' ) {
        ch_filtered_reads = FILTLONG ( reads.map { meta, reads -> [ meta, [], reads ] } ).reads
        ch_versions = ch_versions.mix( FILTLONG.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( FILTLONG.out.log )
    } else if ( params.longread_filter_tool == 'nanoq' ) {
        ch_filtered_reads = NANOQ ( reads , 'fastq.gz' ).reads
        ch_versions        =  ch_versions.mix( NANOQ.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( NANOQ.out.stats )
    } else {
        ch_filtered_reads = reads
    }

    emit:
    reads    = ch_filtered_reads    // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

