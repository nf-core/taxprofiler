//
// Process long raw reads with porechop
//

include { PORECHOP_PORECHOP          } from '../../modules/nf-core/porechop/porechop/main'

workflow LONGREAD_ADAPTERREMOVAL {
    take:
    reads

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (params.longread_adapterremoval_tool == 'porechop') {
        PORECHOP_PORECHOP ( reads )
        ch_processed_reads = PORECHOP_PORECHOP.out.reads
            .map { meta, reads -> [ meta + [single_end: true], reads ] }
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_PORECHOP.out.log )
    } else {
        ch_processed_reads = reads
    }

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
