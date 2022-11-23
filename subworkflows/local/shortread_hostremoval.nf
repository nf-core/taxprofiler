//
// Remove host reads via alignment and export off-target reads
//

include { BOWTIE2_BUILD             } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN             } from '../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS            } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_VIEW             } from '../../modules/nf-core/samtools/view/main'

workflow SHORTREAD_HOSTREMOVAL {
    take:
    reads     // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta
    index     // /path/to/index

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( !params.shortread_hostremoval_index ) {
        ch_bowtie2_index = BOWTIE2_BUILD ( reference ).index
        ch_versions      = ch_versions.mix( BOWTIE2_BUILD.out.versions )
    } else {
        ch_bowtie2_index = index.first()
    }

    BOWTIE2_ALIGN ( reads, ch_bowtie2_index, true, true)
    ch_versions      = ch_versions.mix( BOWTIE2_ALIGN.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix( BOWTIE2_ALIGN.out.log )

    ch_bowtie2_mapped = BOWTIE2_ALIGN.out.bam
        .map {
            meta, reads ->
                [ meta, reads, [] ]
        }

    SAMTOOLS_VIEW ( ch_bowtie2_mapped, [], [] )
    ch_versions      = ch_versions.mix( SAMTOOLS_VIEW.out.versions.first() )

    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.bam )
    ch_versions      = ch_versions.mix( SAMTOOLS_INDEX.out.versions.first() )

    bam_bai = BOWTIE2_ALIGN.out.bam
        .join(SAMTOOLS_INDEX.out.bai, remainder: true)

    SAMTOOLS_STATS ( bam_bai, reference )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_STATS.out.stats )

    emit:
    stats    = SAMTOOLS_STATS.out.stats
    reads    = BOWTIE2_ALIGN.out.fastq   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions               // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

