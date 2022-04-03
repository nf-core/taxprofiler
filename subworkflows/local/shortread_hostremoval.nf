//
// Remove host reads via alignment and export off-target reads
//

include { BOWTIE2_ALIGN             } from '../../../modules/nf-core/modules/bowtie2/align/main'
include { BOWTIE2_BUILD             } from '../../../modules/nf-core/modules/bowtie2/build/main'
include { SAMTOOLS_VIEW             } from '../../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_FASTQ            } from '../../../modules/nf-core/modules/samtools/fastq/main'
include { SAMTOOLS_FLAGSTAT         } from '../../../modules/nf-core/modules/samtools/flagstat/main'

workflow SHORTREAD_PREPROCESSING {
    take:
    reads     // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( !params.shortread_hostremoval_index ) {
        file( , checkIfExists: true )
        BOWTIE2_BUILD ( reference )
        ch_versions = ch_versions.mix( BOWTIE2_BUILD.out.versions )
    }

    BOWTIE2_ALIGN ( reads, BOWTIE2_BUILD.out.index )
    ch_versions = ch_versions.mix( BOWTIE2_BUILD.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT.out.log )

    SAMTOOLS_FLAGSTAT ( BOWTIE2_ALIGN.out.bam )
    ch_versions = ch_versions.mix( SAMTOOLS_FLAGSTAT.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT.out.flagstat )

    emit:
    reads    = BOWTIE2_ALIGN.out.fastq   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

