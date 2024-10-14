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
            .filter { meta, sample_id, instrument_platform,fastq_1,fastq_2,fasta -> (fastq_1 && fastq_2) && !fasta }.view()
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

    channelToSamplesheet(ch_colnames, ch_list_for_samplesheet, 'downstream_samplesheets', format, format_sep)
}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {

    take:
    ch_processed_reads

    main:
    if ( params.generate_pipeline_samplesheets == 'mag' && params.save_analysis_ready_fastqs ) {
        SAMPLESHEET_MAG(ch_processed_reads)
    }
}

def channelToSamplesheet(ch_header, ch_list_for_samplesheet, outdir_subdir, format, format_sep) {
    // Constructs the header string and then the strings of each row, and
    // finally concatenates for saving. Originally designed by @mahesh-panchal
    ch_header
        .first()
        .map { it.keySet().join(format_sep) }
        .concat(ch_list_for_samplesheet.map { it.values().join(format_sep) })
        .collectFile(
            name: "${params.outdir}/${outdir_subdir}/${params.generate_pipeline_samplesheets}.${format}",
            newLine: true,
            sort: false
        )
}
