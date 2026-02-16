include { NONPAREIL_NONPAREIL        } from '../../modules/nf-core/nonpareil/nonpareil/main'
include { NONPAREIL_CURVE            } from '../../modules/nf-core/nonpareil/curve/main'
include { NONPAREIL_SET              } from '../../modules/nf-core/nonpareil/set/main'
include { NONPAREIL_NONPAREILCURVESR } from '../../modules/nf-core/nonpareil/nonpareilcurvesr/main'

workflow NONPAREIL {
    take:
    reads // [ [ meta ], [ reads ] ]

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_reads_for_nonpareil = reads
        .map { meta, input_reads ->
            def reads_new = meta.single_end ? input_reads : input_reads[0]
            // taxprofiler only accepts gzipped input files,
            // so don't need to account for getBaseName removing all extensions
            def format = reads_new[0].getBaseName().split('\\.').last() in ['fasta', 'fna', 'fa', 'fas'] ? 'fasta' : 'fastq'
            [meta, reads_new, format]
        }
        .multiMap { meta, input_reads, format ->
            reads: [meta, input_reads]
            format: format
        }

    // Calculation
    NONPAREIL_NONPAREIL(ch_reads_for_nonpareil.reads, ch_reads_for_nonpareil.format, params.shortread_redundancyestimation_mode)

    ch_npos_for_nonparielset = NONPAREIL_NONPAREIL.out.npo
        .map { _meta, npo -> [[id: 'all'], npo] }
        .groupTuple()

    // Plotting
    NONPAREIL_CURVE(NONPAREIL_NONPAREIL.out.npo)
    // For static single-curve PNG
    NONPAREIL_SET(ch_npos_for_nonparielset)
    // For static multi-curve PNG
    NONPAREIL_NONPAREILCURVESR(ch_npos_for_nonparielset)
    // For dynamic multi-curve PNG in MultiQC and raw files

    ch_versions = ch_versions.mix(
        NONPAREIL_NONPAREIL.out.versions.first(),
        NONPAREIL_CURVE.out.versions.first(),
        NONPAREIL_SET.out.versions.first(),
        NONPAREIL_NONPAREILCURVESR.out.versions.first(),
    )
    ch_multiqc_files = ch_multiqc_files.mix(NONPAREIL_NONPAREILCURVESR.out.json)

    emit:
    npo              = NONPAREIL_NONPAREIL.out.npo
    curve_pngs       = NONPAREIL_CURVE.out.png
    set_pngs         = NONPAREIL_SET.out.png
    nonpareilcurvesr = NONPAREIL_NONPAREILCURVESR.out.tsv
    versions         = ch_versions
    mqc              = ch_multiqc_files
}
