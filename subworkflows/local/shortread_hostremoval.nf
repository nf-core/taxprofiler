//
// Remove host reads via alignment and export off-target reads
//

include { BOWTIE2_BUILD                             } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                             } from '../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX                            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                            } from '../../modules/nf-core/samtools/stats/main'
include { HOSTILE_FETCH as HOSTILE_FETCH_SHORTREADS } from '../../modules/nf-core/hostile/fetch/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_SHORTREADS } from '../../modules/nf-core/hostile/clean/main'

workflow SHORTREAD_HOSTREMOVAL {
    take:
    reads     // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta
    index     // /path/to/index

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (!params.shortread_hostremoval_index && params.shortread_hostremoval_tool == 'bowtie2') {
        ch_hostremoval_index = BOWTIE2_BUILD([[], reference]).index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }
    else if (!params.shortread_hostremoval_index && !index && params.shortread_hostremoval_tool == 'hostile') {
        HOSTILE_FETCH_SHORTREADS(params.hostremoval_hostile_referencename)
        ch_versions = ch_versions.mix(HOSTILE_FETCH_SHORTREADS.out.versions)
        ch_hostremoval_index = HOSTILE_FETCH_SHORTREADS.out.reference
    }
    else {
        ch_hostremoval_index = index.first()
    }

    if (params.shortread_hostremoval_tool == 'bowtie2') {
        // Map, generate BAM with all reads and unmapped reads in FASTQ for downstream
        BOWTIE2_ALIGN(reads, ch_hostremoval_index, [[], reference], true, true)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log)

        ch_cleaned_reads = BOWTIE2_ALIGN.out.fastq

        // Indexing whole BAM for host removal statistics
        SAMTOOLS_INDEX(BOWTIE2_ALIGN.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

        bam_bai = BOWTIE2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai, remainder: true)

        SAMTOOLS_STATS(bam_bai, [[], reference])
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)
    }
    else if (params.shortread_hostremoval_tool == 'hostile') {
        // HOSTILE specifically needs a name of the reference to either download or
        // find the correct files in the index directory
        ch_hostremoval_index_hostile = ch_hostremoval_index.map { _meta, indexdir -> [params.hostremoval_hostile_referencename, indexdir] }

        HOSTILE_CLEAN_SHORTREADS(reads, ch_hostremoval_index_hostile)
        ch_versions = ch_versions.mix(HOSTILE_CLEAN_SHORTREADS.out.versions)
        ch_cleaned_reads = HOSTILE_CLEAN_SHORTREADS.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(HOSTILE_CLEAN_SHORTREADS.out.json)
    }

    emit:
    reads    = ch_cleaned_reads // channel: [ val(meta), [ reads ] ]
    versions = ch_versions // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
