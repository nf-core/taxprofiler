/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTaxprofiler.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.databases, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input    ) { ch_input     = file(params.input)     } else { exit 1, 'Input samplesheet not specified!' }
if (params.databases) { ch_databases = file(params.databases) } else { exit 1, 'Input database sheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'

include { DB_CHECK            } from '../subworkflows/local/db_check'
include { FASTQ_PREPROCESSING } from '../subworkflows/local/preprocessing'


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TAXPROFILER {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    DB_CHECK (
        ch_databases
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: Run Clip/Merge/Complexity
    //
    if ( params.fastp_clip_merge ) {
        FASTQ_PREPROCESSING ( INPUT_CHECK.out.fastq )
    }

    // MODULE: Cat merge runs of same sample
    ch_processed_for_combine = FASTQ_PREPROCESSING.out.reads
        .dump(tag: "prep_for_combine_grouping")
        .map {
            meta, reads ->
            def meta_new = meta.clone()
            meta_new['run_accession'] = 'combined'
            [ meta_new, reads ]
        }
        .groupTuple ( by: 0 )
        .branch{
            combine: it[1].size() >= 2
            skip: it[1].size() < 2
        }

    CAT_FASTQ ( ch_processed_for_combine.combine )

    // Ready for profiling!
    ch_reads_for_profiling = ch_processed_for_combine.skip
                                .dump(tag: "skip_combine")
                                .mix( CAT_FASTQ.out.reads )
                                .dump(tag: "files_for_profiling")

    // Combine reads with possible databases

    ch_reads_for_profiling.combine(DB_CHECK.out.dbs).dump(tag: "reads_plus_db")

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowTaxprofiler.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    if (params.fastp_clip_merge) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_PREPROCESSING.out.mqc)
    }


    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
