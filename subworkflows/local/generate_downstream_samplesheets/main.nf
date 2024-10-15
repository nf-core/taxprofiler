//
// Subworkflow with functionality specific to the nf-core/mag pipeline
//

workflow SAMPLESHEET_MAG {
    take:
    ch_processed_reads

    main:
    format     = 'csv' // most common format in nf-core
    format_sep = ','


    ch_list_for_samplesheet = ch_processed_reads
            .filter { meta, sample_id, instrument_platform,fastq_1,fastq_2,fasta -> (fastq_1 && fastq_2) && !fasta }
                .map {
                        meta, sample_id, instrument_platform,fastq_1,fastq_2,fasta ->
                            def sample        = meta.id
                            def run           = meta.run_accession  //this should be optional
                            def group         = ""
                            def short_reads_1 = file(params.outdir).toString() + '/' + meta.id + '/' + fastq_1.getName()
                            def short_reads_2 = file(params.outdir).toString() + '/' + meta.id + '/' + fastq_2.getName()
                            def long_reads    = ""
                [sample: sample, run: run, group: group, short_reads_1: short_reads_1, short_reads_2: short_reads_2, long_reads: long_reads]
        }
        .tap { ch_colnames }

    channelToSamplesheet(ch_list_for_samplesheet,"${params.outdir}/downstream_samplesheets/mag", format)
}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {

    take:
    ch_processed_reads

    main:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if ( downstreampipeline_names.contains('mag') && params.save_analysis_ready_fastqs) {
        SAMPLESHEET_MAG(ch_processed_reads)
    }

}

// Constructs the header string and then the strings of each row, and
def channelToSamplesheet(ch_list_for_samplesheet, path, format) {
    format_sep = ["csv":",", "tsv":"\t", "txt":"\t"][format]

    ch_header = ch_list_for_samplesheet

    ch_header
        .first()
        .map{ it.keySet().join(format_sep) }
        .concat( ch_list_for_samplesheet.map{ it.values().join(format_sep) })
        .collectFile(
            name:"${path}.${format}",
            newLine: true,
            sort: false
        )
}
