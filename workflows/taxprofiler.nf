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
def checkPathParamList = [ params.input, params.databases, params.shortread_hostremoval_reference,
                            params.shortread_hostremoval_index, params.multiqc_config
                        ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input    ) { ch_input     = file(params.input)     } else { exit 1, 'Input samplesheet not specified!' }
if (params.databases) { ch_databases = file(params.databases) } else { exit 1, 'Input database sheet not specified!' }
if (params.shortread_clipmerge_mergepairs && params.run_malt ) log.warn "[nf-core/taxprofiler] warning: MALT does not accept uncollapsed paired-reads. Pairs will be profiled as separate files."
if (params.shortread_clipmerge_excludeunmerged && !params.shortread_clipmerge_mergepairs) exit 1, "[nf-core/taxprofiler] error: cannot include unmerged reads when merging not turned on. Please specify --shortread_clipmerge_mergepairs"

// TODO Add check if index but no reference exit 1
if (params.shortread_hostremoval && !params.shortread_hostremoval_reference) { exit 1, "[nf-core/taxprofiler] error: --shortread_hostremoval requested but no --shortread_hostremoval_reference FASTA supplied. Check input." }
if (!params.shortread_hostremoval_reference && params.shortread_hostremoval_reference_index) { exit 1, "[nf-core/taxprofiler] error: --shortread_hostremoval_index provided but no --shortread_hostremoval_reference FASTA supplied. Check input." }

if (params.shortread_hostremoval_reference ) { ch_reference       = file(params.shortread_hostremoval_reference) } else { ch_reference       = [] }
if (params.shortread_hostremoval_index     ) { ch_reference_index = file(params.shortread_hostremoval_index    ) } else { ch_reference_index = [] }

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
include { INPUT_CHECK             } from '../subworkflows/local/input_check'

include { DB_CHECK                      } from '../subworkflows/local/db_check'
include { SHORTREAD_PREPROCESSING       } from '../subworkflows/local/shortread_preprocessing'
include { LONGREAD_PREPROCESSING        } from '../subworkflows/local/longread_preprocessing'
include { SHORTREAD_HOSTREMOVAL         } from '../subworkflows/local/shortread_hostremoval'
include { SHORTREAD_COMPLEXITYFILTERING } from '../subworkflows/local/shortread_complexityfiltering'

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
include { METAPHLAN3                  } from '../modules/nf-core/modules/metaphlan3/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TAXPROFILER {

    ch_versions = Channel.empty()

    /*
        SUBWORKFLOW: Read in samplesheet, validate and stage input files
    */
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    DB_CHECK (
        ch_databases
    )

    /*
        MODULE: Run FastQC
    */
    ch_input_for_fastqc = INPUT_CHECK.out.fastq.mix( INPUT_CHECK.out.nanopore )

    FASTQC (
        ch_input_for_fastqc
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */
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
    
    /*
        SUBWORKFLOW: COMPLEXITY FILTERING
    */

    if ( params.shortread_complexityfilter ) {
        ch_shortreads_filtered = SHORTREAD_COMPLEXITYFILTERING ( ch_shortreads_preprocessed ).reads
    } else {
        ch_shortreads_filtered = ch_shortreads_preprocessed
    }

    /*
        SUBWORKFLOW: HOST REMOVAL
    */

    if ( params.shortread_hostremoval ) {
        ch_shortreads_hostremoved = SHORTREAD_HOSTREMOVAL ( ch_shortreads_filtered, ch_reference, ch_reference_index ).reads
        ch_versions = ch_versions.mix(SHORTREAD_HOSTREMOVAL.out.versions.first())
    } else {
        ch_shortreads_hostremoved = ch_shortreads_filtered
    }

    /*
        COMBINE READS WITH POSSIBLE DATABASES
    */

    // e.g. output [DUMP: reads_plus_db] [['id':'2612', 'run_accession':'combined', 'instrument_platform':'ILLUMINA', 'single_end':1], <reads_path>/2612.merged.fastq.gz, ['tool':'malt', 'db_name':'mal95', 'db_params':'"-id 90"'], <db_path>/malt90]
    ch_input_for_profiling = ch_shortreads_hostremoved
            .mix( ch_longreads_preprocessed )
            .combine(DB_CHECK.out.dbs)
            .branch {
                malt:    it[2]['tool'] == 'malt'
                kraken2: it[2]['tool'] == 'kraken2'
                metaphlan3: it[2]['tool'] == 'metaphlan3'
                unknown: true
            }

    /*
        PREPARE PROFILER INPUT CHANNELS
    */

    // We groupTuple to have all samples in one channel for MALT as database
    // loading takes a long time, so we only want to run it once per database
    // TODO document somewhere we only accept illumina short reads for MALT?
    ch_input_for_malt =  ch_input_for_profiling.malt
                            .filter { it[0]['instrument_platform'] == 'ILLUMINA' }
                            .map {
                                it ->
                                    def temp_meta =  [ id: it[2]['db_name']]  + it[2]
                                    def db = it[3]
                                    [ temp_meta, it[1], db ]
                            }
                            .groupTuple(by: [0,2])
                            .multiMap {
                                it ->
                                    reads: [ it[0], it[1].flatten() ]
                                    db: it[2]
                            }

    // We can run Kraken2 one-by-one sample-wise
    ch_input_for_kraken2 =  ch_input_for_profiling.kraken2
                            .multiMap {
                                it ->
                                    reads: [ it[0] + it[2], it[1] ]
                                    db: it[3]
                            }

    ch_input_for_metaphlan3 = ch_input_for_profiling.metaphlan3
                            .multiMap {
                                it ->
                                    reads: [it[0] + it[2], it[1]]
                                    db: it[3]
                            }

    /*
        MODULE: RUN PROFILING
    */
    if ( params.run_malt ) {
        MALT_RUN ( ch_input_for_malt.reads, params.malt_mode, ch_input_for_malt.db )
    }

    if ( params.run_kraken2 ) {
        KRAKEN2_KRAKEN2 ( ch_input_for_kraken2.reads, ch_input_for_kraken2.db  )
    }

    if ( params.run_metaphlan3 ) {
        METAPHLAN3 ( ch_input_for_metaphlan3.reads, ch_input_for_metaphlan3.db )
    }

    /*
        MODULE: MultiQC
    */

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    workflow_summary    = WorkflowTaxprofiler.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    if (params.shortread_clipmerge) {
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
        ch_versions = ch_versions.mix( SHORTREAD_PREPROCESSING.out.versions )
    }

    if (params.longread_clip) {
        ch_multiqc_files = ch_multiqc_files.mix( LONGREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
        ch_versions = ch_versions.mix( LONGREAD_PREPROCESSING.out.versions )
    }

    if (params.shortread_complexityfilter){
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_COMPLEXITYFILTERING.out.mqc.collect{it[1]}.ifEmpty([]) )
        ch_versions = ch_versions.mix( SHORTREAD_COMPLEXITYFILTERING.out.versions )
    }

    if (params.shortread_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_HOSTREMOVAL.out.mqc.collect{it[1]}.ifEmpty([]))
    }

    if (params.run_kraken2) {
        ch_multiqc_files = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.txt.collect{it[1]}.ifEmpty([])  )
        ch_versions = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
    }

    if (params.run_malt) {
        ch_multiqc_files = ch_multiqc_files.mix( MALT_RUN.out.log.collect{it[1]}.ifEmpty([])  )
        ch_versions = ch_versions.mix( MALT_RUN.out.versions.first() )
    }

    // TODO Versions for Karken/MALT not report?
    // TODO create multiQC module for metaphlan
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
