//
// Remove host reads via alignment and export off-target reads
//

include { BOWTIE2_BUILD             } from '../../modules/nf-core/modules/bowtie2/build/main'
include { BOWTIE2_ALIGN             } from '../../modules/nf-core/modules/bowtie2/align/main'

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

    BOWTIE2_ALIGN ( reads, ch_bowtie2_index, true, false )
    ch_versions      = ch_versions.mix( BOWTIE2_ALIGN.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix( BOWTIE2_ALIGN.out.log )

    emit:
    reads    = BOWTIE2_ALIGN.out.fastq   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions               // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

