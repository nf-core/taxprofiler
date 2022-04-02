/*
Process short raw reads with AdapterRemoval
*/

include { ADAPTERREMOVAL as ADAPTERREMOVAL_SINGLE       } from '../../modules/nf-core/modules/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PAIRED       } from '../../modules/nf-core/modules/adapterremoval/main'
include { CAT_FASTQ                                     } from '../../modules/nf-core/modules/cat/fastq/main'
include {
    ENSURE_FASTQ_EXTENSION as ENSURE_FASTQ_EXTENSION1;
    ENSURE_FASTQ_EXTENSION as ENSURE_FASTQ_EXTENSION2;
    ENSURE_FASTQ_EXTENSION as ENSURE_FASTQ_EXTENSION3;
} from '../../modules/local/ensure_fastq_extension'

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

        ENSURE_FASTQ_EXTENSION1(
            Channel.empty().mix(
                ADAPTERREMOVAL_PAIRED.out.collapsed,
                ADAPTERREMOVAL_PAIRED.out.collapsed_truncated,
                ADAPTERREMOVAL_PAIRED.out.singles_truncated,
                ADAPTERREMOVAL_PAIRED.out.pair1_truncated,
                ADAPTERREMOVAL_PAIRED.out.pair2_truncated
            )
            .map { meta, reads ->
                meta.single_end = true
                [meta, reads]
            }
        )

        CAT_FASTQ(
            ENSURE_FASTQ_EXTENSION1.out.reads
                .groupTuple()
        )

        ENSURE_FASTQ_EXTENSION2(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

        ch_adapterremoval_reads_prepped = CAT_FASTQ.out.reads
            .mix(ENSURE_FASTQ_EXTENSION2.out.reads)

    } else if ( params.shortread_clipmerge_mergepairs && params.shortread_clipmerge_excludeunmerged ) {

        ENSURE_FASTQ_EXTENSION1(
            Channel.empty().mix(
                ADAPTERREMOVAL_PAIRED.out.collapsed,
                ADAPTERREMOVAL_PAIRED.out.collapsed_truncated
            )
            .map { meta, reads ->
                meta.single_end = true
                [meta, reads]
            }
        )

        CAT_FASTQ(
            ENSURE_FASTQ_EXTENSION1.out.reads
                .groupTuple()
        )

        ENSURE_FASTQ_EXTENSION2(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

        ch_adapterremoval_reads_prepped = CAT_FASTQ.out.reads
            .mix(ENSURE_FASTQ_EXTENSION2.out.reads)

    } else {

        ENSURE_FASTQ_EXTENSION1(
            ADAPTERREMOVAL_PAIRED.out.pair1_truncated
            .map { meta, reads ->
                meta.single_end = true
                [meta, reads]
            }
        )

        ENSURE_FASTQ_EXTENSION2(
            ADAPTERREMOVAL_PAIRED.out.pair2_truncated
            .map { meta, reads ->
                meta.single_end = true
                [meta, reads]
            }
        )

        ENSURE_FASTQ_EXTENSION3(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

        ch_adapterremoval_reads_prepped = ENSURE_FASTQ_EXTENSION1.out.reads
            .join(ENSURE_FASTQ_EXTENSION2.out.reads)
            .groupTuple()
            .map { meta, pair1, pair2 ->
                meta.single_end = false
                [ meta, [ pair1, pair2 ].flatten() ]
            }
            .mix(ENSURE_FASTQ_EXTENSION3.out.reads)

    }

    ch_versions = ch_versions.mix( ADAPTERREMOVAL_SINGLE.out.versions.first() )
    ch_versions = ch_versions.mix( ADAPTERREMOVAL_PAIRED.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix(
        ADAPTERREMOVAL_PAIRED.out.log.collect{it[1]},
        ADAPTERREMOVAL_SINGLE.out.log.collect{it[1]}
    )

    emit:
    reads    = ch_adapterremoval_reads_prepped  // channel: [ val(meta), [ reads ] ]
    versions = ch_versions  // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

def ensureFastQExtension(row) {
    def (meta, read) = row
    def filename = file(read.parent).resolve("${read.baseName}.fastq.gz")

    read.renameTo(filename.toString())
    return [meta, read]
}
