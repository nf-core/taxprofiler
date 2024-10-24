//
// Subworkflow with functionality specific to the nf-core/mag pipeline
//

workflow SAMPLESHEET_MAG {
    take:
    ch_processed_reads

    main:
    format = 'csv'

    ch_list_for_samplesheet = ch_processed_reads
        .dump()
        .map { meta, reads ->
            def sample = meta.id
            def run = params.perform_runmerging ? '' : meta.run_accession
            def group = "0"
            //this should be optional
            def short_reads_1 = meta.is_fasta ? "" : file(params.outdir).toString() + '/analysis_ready_fastqs/' + reads[0].getName()
            def short_reads_2 = meta.is_fasta || meta.single_end ? "" : file(params.outdir).toString() + '/analysis_ready_fastqs/' + reads[1].getName()
            def long_reads = meta.is_fasta ? file(params.outdir).toString() + '/analysis_ready_fastqs/' + reads[0].getName() : ""

            [sample: sample, run: run, group: group, short_reads_1: short_reads_1, short_reads_2: short_reads_2, long_reads: long_reads]
        }
        .tap { ch_list_for_samplesheet_all }
        .filter { it.short_reads_1 != "" }
        .branch {
            se: it.short_reads_2 == ""
            pe: it.short_reads_2 != ""
            unknown: true
        }

    // Throw a warning that only long reads are not supported yet by MAG
    ch_list_for_samplesheet_all
        .filter { it.long_reads != "" && it.short_reads_1 == "" }
        .collect { log.warn("[nf-core/taxprofiler] WARNING: Standalone long reads are not yet supported by the nf-core/mag pipeline and will not be in present in `mag-*.csv`. Sample: ${it.sample}") }

    channelToSamplesheet(ch_list_for_samplesheet.pe, "${params.outdir}/downstream_samplesheets/mag-pe", format)
    channelToSamplesheet(ch_list_for_samplesheet.se, "${params.outdir}/downstream_samplesheets/mag-se", format)
}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_processed_reads

    main:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if (downstreampipeline_names.contains('mag') && params.save_analysis_ready_fastqs) {
        SAMPLESHEET_MAG(ch_processed_reads)
    }
}

// Constructs the header string and then the strings of each row, and
def channelToSamplesheet(ch_list_for_samplesheet, path, format) {
    def format_sep = [csv: ",", tsv: "\t", txt: "\t"][format]

    def ch_header = ch_list_for_samplesheet

    ch_header
        .first()
        .map { it.keySet().join(format_sep) }
        .concat(ch_list_for_samplesheet.map { it.values().join(format_sep) })
        .collectFile(
            name: "${path}.${format}",
            newLine: true,
            sort: false
        )
}
