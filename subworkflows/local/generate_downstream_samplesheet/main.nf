//
// Subworkflow with functionality specific to the nf-core/createtaxdb pipeline
//

workflow SAMPLESHEET_DIFFERENTIALABUNDANCE {
    take:
    ch_taxpasta

    main:
    format_sep = '\t'

    ch_taxpasta.map { it ->
        def tool_name = it[0]['tool']
        def id = it[0]['id']
        def file_path = it[1]
        def samplesheet_name = file(file_path).getName()

        ch_list_for_samplesheet = Channel
            .fromPath(file_path)
            .splitCsv(sep: format_sep)
            .map { row -> row.drop(1) }
            .flatten()

        ch_colnames = Channel.of('sample')

        channelToSamplesheet(ch_colnames, ch_list_for_samplesheet, 'downstream_samplesheets/differentialabundance', samplesheet_name )
    }
}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_taxpasta

    mai:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if ( downstreampipeline_names.contains('differentialabundance')) {
        SAMPLESHEET_TAXPROFILER(ch_databases)
    }
}

def channelToSamplesheet(ch_header, ch_list_for_samplesheet, outdir_subdir, samplesheet_name) {
    // Constructs the header string and then the strings of each row, and
    // finally concatenates for saving. Originally designed by @mahesh-panchal
    ch_header
        .concat(ch_list_for_samplesheet)
        .collectFile(
            name: "${params.outdir}/${outdir_subdir}/${samplesheet_name}",
            newLine: true,
            sort: false
        )
}
