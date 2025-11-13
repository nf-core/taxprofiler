include { CAT_FASTQ as MERGE_RUNS } from '../../../modules/nf-core/cat/fastq/main'

//
// Subworkflow to create samplesheets for downstream pipelines
//

workflow SAMPLESHEET_METAVAL {
    take:
    ch_shortreads_filtered
    ch_longreads_preprocessed
    ch_profiles
    ch_classifications
    ch_taxpasta

    main:
    format = 'csv' // most common format in nf-core

    // Process reads with run merging
    ch_processed_reads = ch_shortreads_filtered
        .mix(ch_longreads_preprocessed)
        .map { meta, reads -> [meta - meta.subMap('run_accession'), reads] }
        .groupTuple()
        .map { meta, reads -> [meta, reads.flatten()] }
        .branch { meta, reads ->
            cat: (meta.single_end && reads.size() > 1) || (!meta.single_end && reads.size() > 2)
            skip: true
        }

    MERGE_RUNS(ch_processed_reads.cat)
    ch_processed_reads_runmerged = MERGE_RUNS.out.reads
        .mix(ch_processed_reads.skip)
        .map { meta, reads -> [meta, [reads].flatten(), 1, [:]] } // Added empty db_info map

    // Process profiles and taxpasta with database info
    ch_profiles_taxpasta = ch_profiles
        .filter { meta, profiles -> meta.tool in ['kraken2', 'centrifuge', 'diamond'] }
        .map { meta, profiles ->
            def standard_meta = createStandardMeta(meta)
            def db_info = [:]
            db_info["${meta.tool}_db"] = meta.db_name
            [meta.tool, standard_meta, profiles, meta.db_name]
        }
        .combine(
            ch_taxpasta
                .filter { meta, file -> meta.tool in ['kraken2', 'centrifuge', 'diamond'] }
                .map { meta, taxpasta -> [meta.tool, taxpasta] },
            by: [0]
        )
        .map { key, meta, profiles, db_name, taxpasta ->
            def db_info = [:]
            db_info["${key}_db"] = db_name
            [meta, [profiles, taxpasta], 2, db_info]
        }

    // Process classifications with database info
    ch_classifications_dbinfo = ch_classifications
        .filter { meta, file -> meta.tool in ['kraken2', 'centrifuge'] }
        .map { meta, classifications ->
            def standardMeta = createStandardMeta(meta)
            def db_info = [:]
            db_info["${meta.tool}_db"] = meta.db_name
            [standardMeta, classifications, 3, db_info]
        }

    // Combine all channels and create metaval input
    ch_metaval_input = ch_processed_reads_runmerged
        .mix(ch_profiles_taxpasta, ch_classifications_dbinfo)
        .groupTuple()
        .map { meta, files_list, priorities, db_info_list ->
            // Sort by priority and extract files
            def sorted_data = [files_list, priorities, db_info_list].transpose()
                .sort { it[1] }
            def sorted_files = sorted_data.collect { it[0] }
            // Merge all database info maps
            def merged_db_info = [:]
            db_info_list.each { db_map ->
                merged_db_info += db_map
            }
            [meta, sorted_files, merged_db_info]
        }
        .map { meta, files, db_info ->
            def all_files = files.flatten()
            createMetavalStructure(meta, all_files, db_info)
        }

    // Make your samplesheet channel construct here depending on your downstream
    ch_list_for_samplesheet = ch_metaval_input
        .map { it ->
            //single_end, sample, instrument_platform, fastq_1, fastq_2, kraken2_report, kraken2_results, kraken2_taxpasta, kraken2_db_name,
            //centrifuge_report, centrifuge_result, centrifuge_taxpasta, centrifuge_db_name, diamond, diamond_pasta, diamond_db_name
            def sample              = it[1]
            def instrument_platform = it[2]

            // Define fastq_1 based on platform and merged status using ternary operators
            //def fastq_1             = it[2] == "OXFORD_NANOPORE" ?
            //    (it[3].getName().contains("merged") ?
            //        file(params.outdir).toString() + '/filtered_reads_merged/' + it[3].getName() :
            //        file(params.outdir).toString() + '/' + params.longread_filter_tool + '/' + it[3].getName()) :
            //    (it[3].getName().contains("merged") ?
            //        file(params.outdir).toString() + '/filtered_reads_merged/' + it[3].getName() :
            //        file(params.outdir).toString() + '/bbduk/' + it[3].getName())
//
            //// Fix: Check if fastq_2 is empty list before calling getName()
            //def fastq_2 = (!it[0] && it[4] && !(it[4] instanceof List)) ?
            //    (it[2] == "OXFORD_NANOPORE" ?
            //        (it[4].getName().contains("merged") ?
            //            file(params.outdir).toString() + '/filtered_reads_merged/' + it[4].getName() :
            //            file(params.outdir).toString() + '/' + params.longread_filter_tool + '/' + it[4].getName()) :
            //        (it[4].getName().contains("merged") ?
            //            file(params.outdir).toString() + '/filtered_reads_merged/' + it[4].getName() :
            //            file(params.outdir).toString() + '/bbduk/' + it[4].getName())) : ""

            // Determine fastq_1 path
            def fastq1_name = it[3].getName()

            if (it[2] == "OXFORD_NANOPORE") {
                // Long reads
                if (fastq1_name.contains("merged")) {
                    fastq_1 = "${file(params.outdir)}/filtered_reads_merged/${fastq1_name}"
                } else {
                    fastq_1 = "${file(params.outdir)}/${params.longread_filter_tool}/${fastq1_name}"
                }
            } else {
                // Short reads
                if (fastq1_name.contains("merged")) {
                    fastq_1 = "${file(params.outdir)}/filtered_reads_merged/${fastq1_name}"
                } else {
                    fastq_1 = "${file(params.outdir)}/bbduk/${fastq1_name}"
                }
            }

            // Determine fastq_2 path
            if (!it[0] && it[4] && !(it[4] instanceof List)) {
                def fastq2_name = it[4].getName()

                if (it[2] == "OXFORD_NANOPORE") {
                    // Long reads
                    if (fastq2_name.contains("merged")) {
                        fastq_2 = "${file(params.outdir)}/filtered_reads_merged/${fastq2_name}"
                    } else {
                        fastq_2 = "${file(params.outdir)}/${params.longread_filter_tool}/${fastq2_name}"
                    }
                } else {
                    // Short reads
                    if (fastq2_name.contains("merged")) {
                        fastq_2 = "${file(params.outdir)}/filtered_reads_merged/${fastq2_name}"
                    } else {
                        fastq_2 = "${file(params.outdir)}/bbduk/${fastq2_name}"
                    }
                }
            } else {
                fastq_2 = ""
            }

            // Fix: Check if kraken2 files exist before calling getName()
            def kraken2_report      = (it[5] && !(it[5] instanceof List)) ?
                file(params.outdir).toString() + '/kraken2/' + it[8] + '/' + it[5].getName() : ""
            def kraken2_result      = (it[6] && !(it[6] instanceof List)) ?
                file(params.outdir).toString() + '/kraken2/' + it[8] + '/' + it[6].getName() : ""
            def kraken2_taxpasta    = (it[7] && !(it[7] instanceof List)) ?
                file(params.outdir).toString() + '/taxpasta/' + it[7].getName() : ""

            // Fix: Check if centrifuge files exist before calling getName()
            def centrifuge_report   = (it[9] && !(it[9] instanceof List)) ?
                file(params.outdir).toString() + '/centrifuge/' + it[12] + '/' + it[9].getName() : ""
            def centrifuge_result   = (it[10] && !(it[10] instanceof List)) ?
                file(params.outdir).toString() + '/centrifuge/' + it[12] + '/' + it[10].getName() : ""
            def centrifuge_taxpasta = (it[11] && !(it[11] instanceof List)) ?
                file(params.outdir).toString() + '/taxpasta/' + it[11].getName() : ""

            // Fix: Check if diamond files exist before calling getName()
            def diamond = (it[13] && !(it[13] instanceof List)) ?
                file(params.outdir).toString() + '/diamond/' + it[15] + '/' + it[13].getName() : ""
            def diamond_pasta = (it[14] && !(it[14] instanceof List)) ?
                file(params.outdir).toString() + '/taxpasta/' + it[14].getName() : ""

            [ sample: sample, instrument_platform:instrument_platform, fastq_1:fastq_1, fastq_2:fastq_2, kraken2_report: kraken2_report,
            kraken2_result: kraken2_result, kraken2_taxpasta: kraken2_taxpasta, centrifuge_report: centrifuge_report, centrifuge_result: centrifuge_result,
            centrifuge_taxpasta: centrifuge_taxpasta, diamond: diamond, diamond_pasta: diamond_pasta ]
        }

    channelToSamplesheet(ch_list_for_samplesheet,"${params.outdir}/downstream_samplesheets/metaval", format)
}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_shortreads_filtered
    ch_longreads_preprocessed
    ch_profiles
    ch_classifications
    ch_taxpasta

    main:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if ( downstreampipeline_names.contains('metaval')) {
        SAMPLESHEET_METAVAL( ch_shortreads_filtered, ch_longreads_preprocessed, ch_profiles, ch_classifications, ch_taxpasta)
    }
}

// Helper functions

// Create standardized metadata
def createStandardMeta(meta) {
    return [
        id: meta.id,
        instrument_platform: meta.instrument_platform,
        single_end: meta.single_end,
        is_fasta: meta.is_fasta,
        type: meta.type
    ]
}

// Find files by pattern
def findFileByPattern(files, pattern) {
    return files.find { it.toString().contains(pattern) } ?: []
}

// Create metaval structure with database names
def createMetavalStructure(meta, all_files, db_info) {
    return [
        meta.single_end, // 0 - single_end
        meta.id,  // 1 - sample
        meta.instrument_platform,  // 2 - instrument_platform
        all_files[0],  // 3 - fastq_1
        meta.single_end ? [] : (all_files.size() > 1 ? all_files[1] : []),  // 4 - fastq_2
        findFileByPattern(all_files, 'kraken2.report'),  // 5 - kraken2_report
        findFileByPattern(all_files, 'kraken2.classifiedreads'),  // 6 - kraken2_results
        all_files.find { it.toString().contains('kraken2') && it.toString().contains('.tsv') } ?: [], // 7 - kraken2_taxpasta
        db_info.kraken2_db ?: "",  // 8 - kraken2_db_name
        findFileByPattern(all_files, 'centrifuge.txt'),  // 9 - centrifuge_report
        findFileByPattern(all_files, 'centrifuge.results'),  // 10 - centrifuge_result
        all_files.find { it.toString().contains('centrifuge') && it.toString().contains('.tsv') } ?: [],  // 11 - centrifuge_taxpasta
        db_info.centrifuge_db ?: "",  // 12 - centrifuge_db_name
        findFileByPattern(all_files, 'diamond.tsv'),  // 13 - diamond
        all_files.find { it.toString().contains('diamond') && it.toString().contains('.tsv') && !it.toString().contains('diamond.tsv') } ?: [], // 14 - diamond_pasta
        db_info.diamond_db ?: ""  // 15 - diamond_db_name
    ]
}

// Constructs the header string and then the strings of each row
def channelToSamplesheet(ch_list_for_samplesheet, path, format) {
    def format_sep = ["csv":",", "tsv":"\t", "txt":"\t"][format]

    def ch_header = ch_list_for_samplesheet

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
