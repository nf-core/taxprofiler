//
// Remove host reads via alignment and export off-target reads
//

include { BOWTIE2_BUILD                             } from '../../../modules/nf-core/bowtie2/build'
include { BOWTIE2_ALIGN                             } from '../../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_INDEX                            } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_STATS                            } from '../../../modules/nf-core/samtools/stats'
include { HOSTILE_FETCH as HOSTILE_FETCH_SHORTREADS } from '../../../modules/nf-core/hostile/fetch'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_SHORTREADS } from '../../../modules/nf-core/hostile/clean'

workflow SHORTREAD_HOSTREMOVAL {
    take:
    ch_reads // [ [ meta ], [ reads ] ]
    ch_reference // /path/to/fasta
    ch_index // /path/to/index

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    if (reference && !index && params.shortread_hostremoval_tool == 'bowtie2') {
        ch_hostremoval_index = BOWTIE2_BUILD([[], reference]).index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }
    else if (!index && params.shortread_hostremoval_tool == 'hostile') {
        HOSTILE_FETCH_SHORTREADS(params.hostremoval_hostile_referencename)
        ch_versions = ch_versions.mix(HOSTILE_FETCH_SHORTREADS.out.versions)
        ch_hostremoval_index = HOSTILE_FETCH_SHORTREADS.out.reference
    }
    else {
        ch_hostremoval_index = index
    }

    if (params.shortread_hostremoval_tool == 'bowtie2') {
        // Map, generate BAM with all reads and unmapped reads in FASTQ for downstream
        BOWTIE2_ALIGN(reads, ch_hostremoval_index, [[], reference], true, true)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log)

        ch_cleaned_reads = BOWTIE2_ALIGN.out.fastq

        // Indexing whole BAM for host removal statistics
        SAMTOOLS_INDEX(BOWTIE2_ALIGN.out.bam)
        ch_bam_bai = BOWTIE2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.index, remainder: true)

        SAMTOOLS_STATS(ch_bam_bai, [[], reference, []])
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
