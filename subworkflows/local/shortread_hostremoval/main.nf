//
// Remove host reads via alignment and export off-target reads
//

include { BOWTIE2_BUILD  } from '../../../modules/nf-core/bowtie2/build'
include { BOWTIE2_ALIGN  } from '../../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats'

workflow SHORTREAD_HOSTREMOVAL {
    take:
    ch_reads // [ [ meta ], [ reads ] ]
    ch_reference // /path/to/fasta
    ch_index // /path/to/index

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    if (!params.shortread_hostremoval_index) {
        ch_bowtie2_index = BOWTIE2_BUILD([[], ch_reference]).index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }
    else {
        ch_bowtie2_index = ch_index
    }

    // Map, generate BAM with all reads and unmapped reads in FASTQ for downstream
    BOWTIE2_ALIGN(ch_reads, ch_bowtie2_index, [[], ch_reference], true, true)
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log)

    // Indexing whole BAM for host removal statistics
    SAMTOOLS_INDEX(BOWTIE2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    bam_bai = BOWTIE2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai, remainder: true)

    SAMTOOLS_STATS(bam_bai, [[], ch_reference])
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)

    emit:
    stats    = SAMTOOLS_STATS.out.stats
    reads    = BOWTIE2_ALIGN.out.fastq // channel: [ val(meta), [ reads ] ]
    versions = ch_versions // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
