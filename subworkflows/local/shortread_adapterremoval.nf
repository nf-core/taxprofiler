/*
Process short raw reads with AdapterRemoval
*/

include { ADAPTERREMOVAL as ADAPTERREMOVAL_SINGLE       } from '../../modules/nf-core/modules/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PAIRED       } from '../../modules/nf-core/modules/adapterremoval/main'
include { CAT_FASTQ                                     } from '../../modules/nf-core/modules/cat/fastq/main'

workflow SHORTREAD_ADAPTERREMOVAL {

    take:
    reads // [[meta], [reads]]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files      = Channel.empty()

    ch_input_for_adapterremoval = reads
                                    .branch{
                                        single: it[0].single_end
                                        paired: !it[0].single_end
                                    }

    ADAPTERREMOVAL_SINGLE ( ch_input_for_adapterremoval.single, [] )
    ADAPTERREMOVAL_PAIRED ( ch_input_for_adapterremoval.paired, [] )

    /*
     * Due to the ~slightly~ very ugly output implementation of the current AdapterRemoval2 version, each file
     * has to be exported in a separate channel and we must manually recombine when necessary.
     */

    if ( params.shortread_clipmerge_mergepairs && !params.shortread_clipmerge_excludeunmerged ) {

        ch_concat_fastq = Channel.empty()
            .mix(
                ADAPTERREMOVAL_PAIRED.out.collapsed,
                ADAPTERREMOVAL_PAIRED.out.collapsed_truncated,
                ADAPTERREMOVAL_PAIRED.out.singles_truncated,
                ADAPTERREMOVAL_PAIRED.out.paired_truncated
            )
            .groupTuple()
            .map { [it.head(), it.tail().flatten()] }

        CAT_FASTQ(ch_concat_fastq)

        ch_adapterremoval_reads_prepped = CAT_FASTQ.out.reads
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

    } else if ( params.shortread_clipmerge_mergepairs && params.shortread_clipmerge_excludeunmerged ) {

        ch_concat_fastq = Channel.empty()
            .mix(
                ADAPTERREMOVAL_PAIRED.out.collapsed,
                ADAPTERREMOVAL_PAIRED.out.collapsed_truncated
            )
            .map { meta, reads ->
                meta.single_end = true
                [meta, reads]
            }
            .groupTuple()

        CAT_FASTQ(ch_concat_fastq)

        ch_adapterremoval_reads_prepped = CAT_FASTQ.out.reads
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

    } else {

        ch_adapterremoval_reads_prepped = ADAPTERREMOVAL_PAIRED.out.paired_truncated
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

    }

    ch_versions = ch_versions.mix( ADAPTERREMOVAL_SINGLE.out.versions.first() )
    ch_versions = ch_versions.mix( ADAPTERREMOVAL_PAIRED.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix(
        ADAPTERREMOVAL_PAIRED.out.settings.collect{it[1]},
        ADAPTERREMOVAL_SINGLE.out.settings.collect{it[1]}
    )

    emit:
    reads    = ch_adapterremoval_reads_prepped  // channel: [ val(meta), [ reads ] ]
    versions = ch_versions  // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

