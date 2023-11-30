include { NONPAREIL_NONPAREIL } from '../../modules/nf-core/nonpareil/nonpareil/main'
include { NONPAREIL_CURVE     } from '../../modules/nf-core/nonpareil/curve/main'
include { NONPAREIL_SET       } from '../../modules/nf-core/nonpareil/set/main'

// Custom Functions

/*

*/
def extractNonpareilExtensionFromArrays(ch_input) {

return ch_profile
.map { meta, profile -> [meta.db_name, meta, profile] }
    .combine(ch_database, by: 0)
    .multiMap {
        key, meta, profile, db_meta, db ->
            profile: [meta, profile]
            db: db
    }
}

workflow NONPAREIL {
    take:
    reads     // [ [ meta ], [ reads ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    ch_reads_for_nonpareil = reads
                                .map {
                                    meta, reads ->
                                        def reads_new = meta.single_end ? reads : reads[0]
                                        def format = reads_new[0].getBaseName().split('\\.').last() in ['fasta', 'fna', 'fa', 'fas'] ? 'fasta' : 'fastq'

                                    [meta, reads_new, format]
                                }
                                .multiMap {
                                    meta, reads, format ->
                                        reads: [meta, reads]
                                        format: format
                                }

    NONPAREIL_NONPAREIL( ch_reads_for_nonpareil.reads, ch_reads_for_nonpareil.format, params.shortread_redundancyestimation_mode)
    NONPAREIL_CURVE ( NONPAREIL_NONPAREIL.out.npo )
    NONPAREIL_SET ( NONPAREIL_NONPAREIL.out.npo.map {meta, npo -> [[id: 'all'], npo] }.groupTuple() )

    ch_versions = ch_versions.mix( NONPAREIL_NONPAREIL.out.versions.first(), NONPAREIL_CURVE.out.versions.first(), NONPAREIL_SET.out.versions.first() )

    ch_multiqc_files = ch_multiqc_files.mix( NONPAREIL_SET.out.png )

    emit:
    npo        = NONPAREIL_NONPAREIL.out.npo
    curve_pngs = NONPAREIL_CURVE.out.png
    set_pngs   = NONPAREIL_CURVE.out.png
    versions   = ch_versions
    mqc        = ch_multiqc_files
}
