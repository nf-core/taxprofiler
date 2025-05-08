//
// Process long raw reads with porechop or porechop_abi
//

include { PORECHOP_PORECHOP          } from '../../modules/nf-core/porechop/porechop/main'
include { PORECHOP_ABI               } from '../../modules/nf-core/porechop/abi/main'

workflow LONGREAD_ADAPTERREMOVAL {
    take:
    reads

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( params.longread_adapterremoval_tool == 'porechop_abi' ) {
        PORECHOP_ABI ( reads )
        ch_processed_reads = PORECHOP_ABI.out.reads
            .map { meta, reads -> [meta + [ single_end: true ], reads ] }
        ch_versions = ch_versions.mix( PORECHOP_ABI.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_ABI.out.log )
    } else if ( params.longread_adapterremoval_tool == 'porechop' ) {
        PORECHOP_PORECHOP ( reads )
        ch_processed_reads = PORECHOP_PORECHOP.out.reads
            .map { meta, reads -> [ meta + [ single_end: true ], reads ] }
        ch_versions = ch_versions.mix( PORECHOP_PORECHOP.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_PORECHOP.out.log )
    } else {
        ch_processed_reads = reads
    }

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
