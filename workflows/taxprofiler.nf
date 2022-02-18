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
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
include { INPUT_CHECK } from '../subworkflows/local/input_check'

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

include { FASTP as FASTP_SINGLE       } from '../modules/nf-core/modules/fastp/main'
include { FASTP as FASTP_PAIRED       } from '../modules/nf-core/modules/fastp/main'
include { FASTQC as FASTQC_POST       } from '../modules/nf-core/modules/fastqc/main'
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
    // TODO give option to clip only and retain pairs
    // TODO give option to retain singletons (probably fastp option likely)
    // TODO move to subworkflow
    if ( params.fastp_clip_merge ) {

        ch_input_for_fastp = INPUT_CHECK.out.fastq
                                .dump(tag: "pre-fastp_branch")
                                .branch{
                                    single: it[0]['single_end'] == true
                                    paired: it[0]['single_end'] == false
                                }

        ch_input_for_fastp.single.dump(tag: "input_fastp_single")
        ch_input_for_fastp.paired.dump(tag: "input_fastp_paired")

        FASTP_SINGLE ( ch_input_for_fastp.single, false, false )
        FASTP_PAIRED ( ch_input_for_fastp.paired, false, true )

        ch_fastp_reads_prepped = FASTP_PAIRED.out.reads_merged
                                    .mix( FASTP_SINGLE.out.reads )
                                    .map {
                                        meta, reads ->
                                        def meta_new = meta.clone()
                                        meta_new['single_end'] = 1
                                        [ meta_new, reads ]
                                    }

        FASTQC_POST ( ch_fastp_reads_prepped )

        ch_versions = ch_versions.mix(FASTP_SINGLE.out.versions.first())
        ch_versions = ch_versions.mix(FASTP_PAIRED.out.versions.first())

        ch_processed_reads = ch_fastp_reads_prepped

    } else {
        ch_processed_reads = INPUT_CHECK.out.fastq
    }


    // MODULE: Cat merge runs of same sample
    ch_processed_for_combine = ch_processed_reads
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
        ch_multiqc_files = ch_multiqc_files.mix(FASTP_SINGLE.out.json.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTP_PAIRED.out.json.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect{it[1]}.ifEmpty([]))
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
