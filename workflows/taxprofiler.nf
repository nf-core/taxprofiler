/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'

include { DB_CHECK            } from '../subworkflows/local/db_check'
include { SHORTREAD_PREPROCESSING } from '../subworkflows/local/shortread_preprocessing'
include { LONGREAD_PREPROCESSING } from '../subworkflows/local/longread_preprocessing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { MALT_RUN                    } from '../modules/nf-core/modules/malt/run/main'
include { KRAKEN2_KRAKEN2             } from '../modules/nf-core/modules/kraken2/kraken2/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    ch_input_for_fastqc = INPUT_CHECK.out.fastq.mix( INPUT_CHECK.out.nanopore ).dump(tag: "input_to_fastq")
    FASTQC (
        ch_input_for_fastqc
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // PERFORM PREPROCESSING
    //
    if ( params.shortread_clipmerge ) {
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING ( INPUT_CHECK.out.fastq ).reads
    } else {
        ch_shortreads_preprocessed = INPUT_CHECK.out.fastq
    }

    if ( params.longread_clip ) {
        ch_longreads_preprocessed = LONGREAD_PREPROCESSING ( INPUT_CHECK.out.nanopore ).reads
                                        .map { it -> [ it[0], [it[1]] ] }
    ch_versions = ch_versions.mix(LONGREAD_PREPROCESSING.out.versions.first())
    } else {
        ch_longreads_preprocessed = INPUT_CHECK.out.nanopore
    }

    //
    // PERFORM SHORT READ RUN MERGING
    // TODO: Check not necessary for long reads too?
    //
    ch_processed_for_combine = ch_shortreads_preprocessed
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

    ch_reads_for_profiling = ch_processed_for_combine.skip
                                .dump(tag: "skip_combine")
                                .mix( CAT_FASTQ.out.reads )
                                .dump(tag: "files_for_profiling")

    //
    // COMBINE READS WITH POSSIBLE DATABASES
    //

    // e.g. output [DUMP: reads_plus_db] [['id':'2612', 'run_accession':'combined', 'instrument_platform':'ILLUMINA', 'single_end':1], <reads_path>/2612.merged.fastq.gz, ['tool':'malt', 'db_name':'mal95', 'db_params':'"-id 90"'], <db_path>/malt90]
    ch_input_for_profiling = ch_reads_for_profiling
            .mix( ch_longreads_preprocessed )
            .combine(DB_CHECK.out.dbs)
            .dump(tag: "reads_plus_db")
            .branch {
                malt:    it[2]['tool'] == 'malt'
                kraken2: it[2]['tool'] == 'kraken2'
                unknown: true
            }

    //
    // PREPARE PROFILER INPUT CHANNELS
    //

    // We groupTuple to have all samples in one channel for MALT as database
    // loading takes a long time, so we only want to run it once per database
    ch_input_for_malt =  ch_input_for_profiling.malt
                            .map {
                                it ->
                                    def temp_meta =  [ id: it[2]['db_name']]  + it[2]
                                    def db = it[3]
                                    [ temp_meta, it[1], db ]
                            }
                            .groupTuple(by: [0,2])
                            .dump(tag: "input for malt")
                            .multiMap {
                                it ->
                                    reads: [ it[0], it[1].flatten() ]
                                    db: it[2]
                            }

    // We can run Kraken2 one-by-one sample-wise
    ch_input_for_kraken2 =  ch_input_for_profiling.kraken2
                            .dump(tag: "input for kraken")
                            .multiMap {
                                it ->
                                    reads: [ it[0] + it[2], it[1] ]
                                    db: it[3]
                            }

    //
    // RUN PROFILING
    //
    if ( params.run_malt ) {
        MALT_RUN ( ch_input_for_malt.reads, params.malt_mode, ch_input_for_malt.db )
    }

    if ( params.run_kraken2 ) {
        KRAKEN2_KRAKEN2 ( ch_input_for_kraken2.reads, ch_input_for_kraken2.db  )
    }

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

    if (params.shortread_clipmerge) {
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_PREPROCESSING.out.mqc)
    }
    if (params.longread_clip) {
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_PREPROCESSING.out.mqc)
    }
    if (params.run_kraken2) {
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.txt.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())
    }
    if (params.run_malt) {
        ch_multiqc_files = ch_multiqc_files.mix(MALT_RUN.out.log.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(MALT_RUN.out.versions.first())
    }

    // TODO MALT results overwriting per database?
    // TODO Versions for Karken/MALT not report?
    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
