//
// Remove host reads via alignment and export off-target reads
//

include { MINIMAP2_INDEX             } from '../../modules/nf-core/modules/minimap2/index/main'
include { MINIMAP2_ALIGN             } from '../../modules/nf-core/modules/minimap2/align/main'
include { SAMTOOLS_VIEW              } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_BAM2FQ            } from '../../modules/nf-core/modules/samtools/bam2fq/main'

workflow LONGREAD_HOSTREMOVAL {
    take:
    reads     // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta
    index     // /path/to/index

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( !params.longread_hostremoval_index ) {
        ch_minimap2_index = MINIMAP2_INDEX ( reference ).index
        ch_versions       = ch_versions.mix( MINIMAP2_INDEX.out.versions )
    } else {
        ch_minimap2_index = index
    }

    MINIMAP2_ALIGN ( reads, ch_minimap2_index, true, false, false )
    ch_versions        = ch_versions.mix( MINIMAP2_ALIGN.out.versions.first() )
    ch_minimap2_mapped = MINIMAP2_ALIGN.out.bam
        .map {
            meta, reads ->
                [ meta, reads, [] ]
        }


    SAMTOOLS_VIEW ( ch_minimap2_mapped , [] )
    ch_versions      = ch_versions.mix( SAMTOOLS_VIEW.out.versions.first() )

    SAMTOOLS_BAM2FQ ( SAMTOOLS_VIEW.out.bam, false )
    ch_versions      = ch_versions.mix( SAMTOOLS_BAM2FQ.out.versions.first() )


    emit:
    reads    = SAMTOOLS_BAM2FQ.out.reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions                 // channel: [ versions.yml ]
}

