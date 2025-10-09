workflow SAMPLESHEET_METAVAL {
    take:
    reads_index
    classifications_index
    profiles_index
    taxpasta_index

    main:
    //format = 'csv' // most common format in nf-core

    ch_reads = Channel.fromPath ( reads_index, checkIfExists: true )
        .splitCsv ( sep:',', header: true )
    ch_classifications = Channel.fromPath( classifications_index, checkIfExists: true )
        .splitCsv ( sep: ',' ,header: true )
    ch_profiles = Channel.fromPath( profiles_index, checkIfExists: true )
        .splitCsv( sep: ',', header: true)
    ch_taxpasta = Channel.fromPath( taxpasta_index, checkIfExists: true)
        .splitCsv( sep: ',', header: true)

    ch_reads.dump (tag:"reads")
    ch_classifications.dump(tag:"classification")
    ch_profiles.dump(tag:"profiles")
    ch_taxpasta.dump(tag:"taxpasta")

    emit:
    reads = ch_reads
    classifications = ch_classifications
    profiles = ch_profiles
    taxpasta = ch_taxpasta

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

