//
// Remove host reads via alignment and export off-target reads
//

include { MINIMAP2_INDEX                           } from '../../../modules/nf-core/minimap2/index'
include { MINIMAP2_ALIGN                           } from '../../../modules/nf-core/minimap2/align'
include { SAMTOOLS_VIEW                            } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_FASTQ                           } from '../../../modules/nf-core/samtools/fastq'
include { SAMTOOLS_STATS                           } from '../../../modules/nf-core/samtools/stats'
include { HOSTILE_FETCH as HOSTILE_FETCH_LONGREADS } from '../../../modules/nf-core/hostile/fetch'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_LONGREADS } from '../../../modules/nf-core/hostile/clean'

workflow LONGREAD_HOSTREMOVAL {
    take:
    reads // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta
    index // /path/to/index

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    if (!params.longread_hostremoval_index && params.longread_hostremoval_tool == 'minimap2') {
        ch_hostremoval_index = MINIMAP2_INDEX([[], reference]).index
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }
    else if (!params.longread_hostremoval_index && params.longread_hostremoval_tool == 'hostile') {
        HOSTILE_FETCH_LONGREADS(params.hostremoval_hostile_referencename)
        ch_versions = ch_versions.mix(HOSTILE_FETCH_LONGREADS.out.versions)
        ch_hostremoval_index = HOSTILE_FETCH_LONGREADS.out.reference
    }
    else {
        ch_hostremoval_index = index
    }

    if (params.longread_hostremoval_tool == 'minimap2') {
        MINIMAP2_ALIGN(reads, ch_hostremoval_index, true, 'bai', false, false)
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
        ch_minimap2_mapped = MINIMAP2_ALIGN.out.bam.map { meta, long_reads ->
            [meta, long_reads, []]
        }

        // Join BAI back to BAM for host removal statistics
        bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index)

        SAMTOOLS_STATS(bam_bai, [[], reference, []])
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)

        // Generate unmapped reads FASTQ for downstream taxprofiling
        SAMTOOLS_VIEW(ch_minimap2_mapped, [[], [], []], [[], []], [[], []], 'bai')
        SAMTOOLS_FASTQ(SAMTOOLS_VIEW.out.bam, false)

        ch_cleaned_reads = SAMTOOLS_FASTQ.out.other
    }
    else if (params.longread_hostremoval_tool == 'hostile') {
        // HOSTILE specifically needs a name of the reference to either download or
        // find the correct files in the index directory
        ch_hostremoval_index_hostile = ch_hostremoval_index.map { _meta, indexdir -> [params.hostremoval_hostile_referencename, indexdir] }

        HOSTILE_CLEAN_LONGREADS(reads, ch_hostremoval_index_hostile)
        ch_versions = ch_versions.mix(HOSTILE_CLEAN_LONGREADS.out.versions)
        ch_cleaned_reads = HOSTILE_CLEAN_LONGREADS.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(HOSTILE_CLEAN_LONGREADS.out.json)
    }

    emit:
    reads    = ch_cleaned_reads // channel: [ val(meta), [ reads ] ]
    versions = ch_versions // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
