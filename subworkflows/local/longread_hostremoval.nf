//
// Remove host reads via alignment and export off-target reads
//

include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_STATS } from '../../modules/nf-core/samtools/stats/main'

workflow LONGREAD_HOSTREMOVAL {
    take:
    reads // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta
    index // /path/to/index

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    if (!params.longread_hostremoval_index) {
        ch_minimap2_index = MINIMAP2_INDEX([[], reference]).index
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }
    else {
        ch_minimap2_index = index
    }

    MINIMAP2_ALIGN(reads, ch_minimap2_index, true, 'bai', false, false)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
    ch_minimap2_mapped = MINIMAP2_ALIGN.out.bam.map { meta, long_reads ->
        [meta, long_reads, []]
    }

    // Join BAI back to BAM for host removal statistics
    bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index)

    SAMTOOLS_STATS(bam_bai, [[], reference])
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)

    // Generate unmapped reads FASTQ for downstream taxprofiling
    SAMTOOLS_VIEW(ch_minimap2_mapped, [[], []], [], 'bai')
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    SAMTOOLS_FASTQ(SAMTOOLS_VIEW.out.bam, false)
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    emit:
    stats    = SAMTOOLS_STATS.out.stats //channel: [val(meta), [reads  ] ]
    reads    = SAMTOOLS_FASTQ.out.other // channel: [ val(meta), [ reads ] ]
    versions = ch_versions // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
