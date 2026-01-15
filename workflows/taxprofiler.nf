/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap              } from 'plugin/nf-schema'
include { paramsSummaryMultiqc          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../subworkflows/local/utils_nfcore_taxprofiler_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { SHORTREAD_PREPROCESSING       } from '../subworkflows/local/shortread_preprocessing'
include { NONPAREIL                     } from '../subworkflows/local/nonpareil'
include { LONGREAD_PREPROCESSING        } from '../subworkflows/local/longread_preprocessing'
include { SHORTREAD_HOSTREMOVAL         } from '../subworkflows/local/shortread_hostremoval'
include { LONGREAD_HOSTREMOVAL          } from '../subworkflows/local/longread_hostremoval'
include { SHORTREAD_COMPLEXITYFILTERING } from '../subworkflows/local/shortread_complexityfiltering'
include { PROFILING                     } from '../subworkflows/local/profiling'
include { VISUALIZATION_KRONA           } from '../subworkflows/local/visualization_krona'
include { STANDARDISATION_PROFILES      } from '../subworkflows/local/standardisation_profiles'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR                         } from '../modules/nf-core/untar/main'
include { FALCO                         } from '../modules/nf-core/falco/main'
include { CAT_FASTQ as MERGE_RUNS       } from '../modules/nf-core/cat/fastq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TAXPROFILER {
    take:
    samplesheet // channel: samplesheet read in from --input
    databases // channel: databases from --databases

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Preprocessing auxiliary file input channel preperation
    adapterlist = params.shortread_qc_adapterlist ? file(params.shortread_qc_adapterlist) : []
    if (params.hostremoval_reference) {
        ch_reference = file(params.hostremoval_reference)
    }
    if (params.shortread_hostremoval_index) {
        ch_shortread_reference_index = channel.fromPath(params.shortread_hostremoval_index).map { [[], it] }
    }
    else {
        ch_shortread_reference_index = []
    }
    if (params.longread_hostremoval_index) {
        ch_longread_reference_index = file(params.longread_hostremoval_index)
    }
    else {
        ch_longread_reference_index = []
    }

    // Validate input files and create separate channels for FASTQ, FASTA, and Nanopore data
    ch_input = samplesheet
        .branch { meta, _run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            fastq: meta.single_end || fastq_2
            return [meta + [type: "short"], fastq_2 ? [fastq_1, fastq_2] : [fastq_1]]
            nanopore: instrument_platform == 'OXFORD_NANOPORE' && !meta.is_fasta
            meta.single_end = true
            return [meta + [type: "long"], [fastq_1]]
            pacbio: instrument_platform == 'PACBIO_SMRT' && !meta.is_fasta
            meta.single_end = true
            return [meta + [type: "long"], [fastq_1]]
            fasta_short: meta.is_fasta && instrument_platform == 'ILLUMINA'
            meta.single_end = true
            return [meta + [type: "short"], [fasta]]
            fasta_long: meta.is_fasta && (instrument_platform == 'OXFORD_NANOPORE' || instrument_platform == 'PACBIO_SMRT')
            meta.single_end = true
            return [meta + [type: "long"], [fasta]]
        }

    // Merge ch_input.fastq and ch_input.nanopore into a single channel
    ch_input_for_fastqc = ch_input.fastq.mix(ch_input.nanopore, ch_input.pacbio)

    // Validate and decompress databases
    ch_dbs_for_untar = databases.branch { db_meta, db_path ->
        if (!db_meta.db_type) {
            db_meta = db_meta + [db_type: "short;long"]
        }
        untar: db_path.name.endsWith(".tar.gz")
        skip: true
    }
    // Filter the channel to untar only those databases for tools that are selected to be run by the user.
    // Also, to ensure only untar once per file, group together all databases of one file
    ch_inputdb_untar = ch_dbs_for_untar.untar
        .filter { db_meta, _db_path ->
            params["run_${db_meta.tool}"]
        }
        .groupTuple(by: 1)
        .map { meta, dbfile ->
            def new_meta = ['id': dbfile.baseName] + ['meta': meta]
            [new_meta, dbfile]
        }

    // Untar the databases
    UNTAR(ch_inputdb_untar)
    ch_versions = ch_versions.mix(UNTAR.out.versions.first())

    // Spread out the untarred and shared databases
    ch_outputdb_from_untar = UNTAR.out.untar
        .map { meta, db ->
            [meta.meta, db]
        }
        .transpose(by: 0)

    ch_final_dbs = ch_dbs_for_untar.skip
        .mix(ch_outputdb_from_untar)
        .map { db_meta, db ->
            def corrected_db_params = db_meta.db_params ? [db_params: db_meta.db_params] : [db_params: '']
            [db_meta + corrected_db_params, db]
        }

    /*
        MODULE: Run FastQC
    */


    if (!params.skip_preprocessing_qc) {
        if (params.preprocessing_qc_tool == 'falco') {
            FALCO(ch_input_for_fastqc)
            ch_versions = ch_versions.mix(FALCO.out.versions.first())
        }
        else {
            FASTQC(ch_input_for_fastqc)
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        }
    }

    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */

    if (params.perform_shortread_qc) {
        SHORTREAD_PREPROCESSING(ch_input.fastq, adapterlist)
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING.out.reads
        ch_versions = ch_versions.mix(SHORTREAD_PREPROCESSING.out.versions)
    }
    else {
        ch_shortreads_preprocessed = ch_input.fastq
    }

    if (params.perform_longread_qc) {
        ch_longreads_preprocessed_nanopore = LONGREAD_PREPROCESSING(ch_input.nanopore).reads.map { it -> [it[0], [it[1]]] }
        ch_versions = ch_versions.mix(LONGREAD_PREPROCESSING.out.versions)
    }
    else {
        ch_longreads_preprocessed_nanopore = ch_input.nanopore
    }

    /*
        MODULE: REDUNDANCY ESTIMATION
    */

    if (params.perform_shortread_redundancyestimation) {
        NONPAREIL(ch_shortreads_preprocessed)
        ch_versions = ch_versions.mix(NONPAREIL.out.versions)
    }

    /*
        SUBWORKFLOW: COMPLEXITY FILTERING
    */

    // fastp complexity filtering is activated via modules.conf in shortread_preprocessing
    if (params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp') {
        SHORTREAD_COMPLEXITYFILTERING(ch_shortreads_preprocessed)
        ch_shortreads_filtered = SHORTREAD_COMPLEXITYFILTERING.out.reads
        ch_versions = ch_versions.mix(SHORTREAD_COMPLEXITYFILTERING.out.versions)
    }
    else {
        ch_shortreads_filtered = ch_shortreads_preprocessed
    }

    /*
        SUBWORKFLOW: HOST REMOVAL
    */

    if (params.perform_shortread_hostremoval) {
        SHORTREAD_HOSTREMOVAL(ch_shortreads_filtered, ch_reference, ch_shortread_reference_index)
        ch_shortreads_hostremoved = SHORTREAD_HOSTREMOVAL.out.reads
        ch_versions = ch_versions.mix(SHORTREAD_HOSTREMOVAL.out.versions)
    }
    else {
        ch_shortreads_hostremoved = ch_shortreads_filtered
    }
    ch_longreads_preprocessed = ch_longreads_preprocessed_nanopore.mix(ch_input.pacbio)
    if (params.perform_longread_hostremoval) {
        LONGREAD_HOSTREMOVAL(ch_longreads_preprocessed, ch_reference, ch_longread_reference_index)
        ch_longreads_hostremoved = LONGREAD_HOSTREMOVAL.out.reads
        ch_versions = ch_versions.mix(LONGREAD_HOSTREMOVAL.out.versions)
    }
    else {
        ch_longreads_hostremoved = ch_longreads_preprocessed
    }

    if (params.perform_runmerging) {

        ch_reads_for_cat_branch = ch_shortreads_hostremoved
            .mix(ch_longreads_hostremoved)
            .map { meta, reads ->
                def meta_new = meta - meta.subMap('run_accession')
                [meta_new, reads]
            }
            .groupTuple()
            .map { meta, reads ->
                [meta, reads.flatten()]
            }
            .branch { meta, reads ->
                cat: (meta.single_end && reads.size() > 1) || (!meta.single_end && reads.size() > 2)
                skip: true
            }
        // we can't concatenate files if there is not a second run, so we branch
        // here to separate them out, and mix back in after for efficiency
        MERGE_RUNS(ch_reads_for_cat_branch.cat)
        ch_reads_runmerged = MERGE_RUNS.out.reads
            .mix(ch_reads_for_cat_branch.skip)
            .map { meta, reads ->
                [meta, [reads].flatten()]
            }
            .mix(ch_input.fasta_short, ch_input.fasta_long)

        ch_versions = ch_versions.mix(MERGE_RUNS.out.versions)
    }
    else {
        ch_reads_runmerged = ch_shortreads_hostremoved.mix(ch_longreads_hostremoved, ch_input.fasta_short, ch_input.fasta_long)
    }

    /*
        SUBWORKFLOW: PROFILING
    */

    PROFILING(ch_reads_runmerged, ch_final_dbs)
    ch_versions = ch_versions.mix(PROFILING.out.versions)

    /*
        SUBWORKFLOW: VISUALIZATION_KRONA
    */
    if (params.run_krona) {
        VISUALIZATION_KRONA(PROFILING.out.classifications, PROFILING.out.profiles, ch_final_dbs)
        ch_versions = ch_versions.mix(VISUALIZATION_KRONA.out.versions)
    }

    /*
        SUBWORKFLOW: PROFILING STANDARDISATION
    */
    if (params.run_profile_standardisation) {
        STANDARDISATION_PROFILES(PROFILING.out.classifications, PROFILING.out.profiles, ch_final_dbs, PROFILING.out.motus_version)
        ch_versions = ch_versions.mix(STANDARDISATION_PROFILES.out.versions)
    }

    /*
        MODULE: MultiQC
    */

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'taxprofiler_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    if (!params.skip_preprocessing_qc) {
        if (params.preprocessing_qc_tool == 'falco') {
            // only mix in files actually used by MultiQC
            ch_multiqc_files = ch_multiqc_files.mix(
                FALCO.out.txt.map { _meta, reports -> reports }.flatten().filter { path -> path.name.endsWith('_data.txt') }.ifEmpty([])
            )
        }
        else {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))
        }
    }

    if (params.perform_shortread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_PREPROCESSING.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    if (params.perform_longread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_PREPROCESSING.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    if (params.perform_shortread_redundancyestimation) {
        ch_multiqc_files = ch_multiqc_files.mix(NONPAREIL.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    if (params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp') {
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_COMPLEXITYFILTERING.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    if (params.perform_shortread_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_HOSTREMOVAL.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    if (params.perform_longread_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_HOSTREMOVAL.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    ch_multiqc_files = ch_multiqc_files.mix(PROFILING.out.mqc.collect { it[1] }.ifEmpty([]))

    if (params.run_profile_standardisation) {
        ch_multiqc_files = ch_multiqc_files.mix(STANDARDISATION_PROFILES.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
