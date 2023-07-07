/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowTaxprofiler.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.genome, params.databases,
                            params.longread_hostremoval_index,
                            params.hostremoval_reference, params.shortread_hostremoval_index,
                            params.multiqc_config, params.shortread_qc_adapterlist,
                            params.krona_taxonomy_directory,
                            params.taxpasta_taxonomy_dir,
                            params.multiqc_logo, params.multiqc_methods_description
                        ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if ( params.input ) {
    ch_input              = file(params.input, checkIfExists: true)
} else {
    error("Input samplesheet not specified")
}

if (params.databases) { ch_databases = file(params.databases, checkIfExists: true) } else { error('Input database sheet not specified!') }

if (!params.shortread_qc_mergepairs && params.run_malt ) log.warn "[nf-core/taxprofiler] MALT does not accept uncollapsed paired-reads. Pairs will be profiled as separate files."
if (params.shortread_qc_includeunmerged && !params.shortread_qc_mergepairs) error("ERROR: [nf-core/taxprofiler] cannot include unmerged reads when merging is not turned on. Please specify --shortread_qc_mergepairs")

if (params.shortread_complexityfilter_tool == 'fastp' && ( params.perform_shortread_qc == false || params.shortread_qc_tool != 'fastp' ))  error("ERROR: [nf-core/taxprofiler] cannot use fastp complexity filtering if preprocessing not turned on and/or tool is not fastp. Please specify --perform_shortread_qc and/or --shortread_qc_tool 'fastp'")

if (params.perform_shortread_hostremoval && !params.hostremoval_reference) { error("ERROR: [nf-core/taxprofiler] --shortread_hostremoval requested but no --hostremoval_reference FASTA supplied. Check input.") }
if (params.perform_shortread_hostremoval && !params.hostremoval_reference && params.shortread_hostremoval_index) { error("ERROR: [nf-core/taxprofiler] --shortread_hostremoval_index provided but no --hostremoval_reference FASTA supplied. Check input.") }
if (params.perform_longread_hostremoval && !params.hostremoval_reference && params.longread_hostremoval_index) { error("ERROR: [nf-core/taxprofiler] --longread_hostremoval_index provided but no --hostremoval_reference FASTA supplied. Check input.") }

if (params.hostremoval_reference           ) { ch_reference = file(params.hostremoval_reference) }
if (params.shortread_hostremoval_index     ) { ch_shortread_reference_index = Channel.fromPath(params.shortread_hostremoval_index).map{[[], it]} } else { ch_shortread_reference_index = [] }
if (params.longread_hostremoval_index      ) { ch_longread_reference_index  = file(params.longread_hostremoval_index     ) } else { ch_longread_reference_index  = [] }

if (params.diamond_save_reads              ) log.warn "[nf-core/taxprofiler] DIAMOND only allows output of a single format. As --diamond_save_reads supplied, only aligned reads in SAM format will be produced, no taxonomic profiles will be available."

if (params.run_malt && params.run_krona && !params.krona_taxonomy_directory) log.warn "[nf-core/taxprofiler] Krona can only be run on MALT output if path to Krona taxonomy database supplied to --krona_taxonomy_directory. Krona will not be executed in this run for MALT."
if (params.run_bracken && !params.run_kraken2) error('ERROR: [nf-core/taxprofiler] You are attempting to run Bracken without running kraken2. This is not possible! Please set --run_kraken2 as well.')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                   } from '../subworkflows/local/input_check'

include { DB_CHECK                      } from '../subworkflows/local/db_check'
include { SHORTREAD_PREPROCESSING       } from '../subworkflows/local/shortread_preprocessing'
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
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { FALCO                       } from '../modules/nf-core/falco/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ as MERGE_RUNS     } from '../modules/nf-core/cat/fastq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TAXPROFILER {

    ch_versions = Channel.empty()
    ch_multiqc_logo= Channel.fromPath("$projectDir/docs/images/nf-core-taxprofiler_logo_custom_light.png")
    adapterlist = params.shortread_qc_adapterlist ? file(params.shortread_qc_adapterlist) : []

    if ( params.shortread_qc_adapterlist ) {
        if ( params.shortread_qc_tool == 'adapterremoval' && !(adapterlist.extension == 'txt') ) error "[nf-core/taxprofiler] ERROR: AdapterRemoval2 adapter list requires a `.txt` format and extension. Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
        if ( params.shortread_qc_tool == 'fastp' && !adapterlist.extension.matches(".*(fa|fasta|fna|fas)") ) error "[nf-core/taxprofiler] ERROR: fastp adapter list requires a `.fasta` format and extension (or fa, fas, fna). Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
    }

    /*
        SUBWORKFLOW: Read in samplesheet, validate and stage input files
    */
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    // Save final FASTA reads if requested, as otherwise no processing occurs on FASTA


    DB_CHECK (
        ch_databases
    )
    ch_versions = ch_versions.mix(DB_CHECK.out.versions)

    /*
        MODULE: Run FastQC
    */
    ch_input_for_fastqc = INPUT_CHECK.out.fastq.mix( INPUT_CHECK.out.nanopore )

    if ( params.preprocessing_qc_tool == 'falco' ) {
        FALCO ( ch_input_for_fastqc )
        ch_versions = ch_versions.mix(FALCO.out.versions.first())
    } else {
        FASTQC ( ch_input_for_fastqc )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }
    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */

    if ( params.perform_shortread_qc ) {
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING ( INPUT_CHECK.out.fastq, adapterlist ).reads
        ch_versions = ch_versions.mix( SHORTREAD_PREPROCESSING.out.versions )
    } else {
        ch_shortreads_preprocessed = INPUT_CHECK.out.fastq
    }

    if ( params.perform_longread_qc ) {
        ch_longreads_preprocessed = LONGREAD_PREPROCESSING ( INPUT_CHECK.out.nanopore ).reads
                                        .map { it -> [ it[0], [it[1]] ] }
        ch_versions = ch_versions.mix( LONGREAD_PREPROCESSING.out.versions )
    } else {
        ch_longreads_preprocessed = INPUT_CHECK.out.nanopore
    }

    /*
        SUBWORKFLOW: COMPLEXITY FILTERING
    */

    // fastp complexity filtering is activated via modules.conf in shortread_preprocessing
    if ( params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp' ) {
        ch_shortreads_filtered = SHORTREAD_COMPLEXITYFILTERING ( ch_shortreads_preprocessed ).reads
        ch_versions = ch_versions.mix( SHORTREAD_COMPLEXITYFILTERING.out.versions )
    } else {
        ch_shortreads_filtered = ch_shortreads_preprocessed
    }

    /*
        SUBWORKFLOW: HOST REMOVAL
    */

    if ( params.perform_shortread_hostremoval ) {
        ch_shortreads_hostremoved = SHORTREAD_HOSTREMOVAL ( ch_shortreads_filtered, ch_reference, ch_shortread_reference_index ).reads
        ch_versions = ch_versions.mix(SHORTREAD_HOSTREMOVAL.out.versions)
    } else {
        ch_shortreads_hostremoved = ch_shortreads_filtered
    }

    if ( params.perform_longread_hostremoval ) {
        ch_longreads_hostremoved = LONGREAD_HOSTREMOVAL ( ch_longreads_preprocessed, ch_reference, ch_longread_reference_index ).reads
        ch_versions = ch_versions.mix(LONGREAD_HOSTREMOVAL.out.versions)
    } else {
        ch_longreads_hostremoved = ch_longreads_preprocessed
    }

    if ( params.perform_runmerging ) {

        ch_reads_for_cat_branch = ch_shortreads_hostremoved
            .mix( ch_longreads_hostremoved )
            .map {
                meta, reads ->
                    def meta_new = meta.clone()
                    meta_new.remove('run_accession')
                    [ meta_new, reads ]
            }
            .groupTuple()
            .map {
                meta, reads ->
                    [ meta, reads.flatten() ]
            }
            .branch {
                meta, reads ->
                // we can't concatenate files if there is not a second run, we branch
                // here to separate them out, and mix back in after for efficiency
                cat: ( meta.single_end && reads.size() > 1 ) || ( !meta.single_end && reads.size() > 2 )
                skip: true
            }

        ch_reads_runmerged = MERGE_RUNS ( ch_reads_for_cat_branch.cat ).reads
            .mix( ch_reads_for_cat_branch.skip )
            .map {
                meta, reads ->
                [ meta, [ reads ].flatten() ]
            }
            .mix( INPUT_CHECK.out.fasta )

        ch_versions = ch_versions.mix(MERGE_RUNS.out.versions)

    } else {
        ch_reads_runmerged = ch_shortreads_hostremoved
            .mix( ch_longreads_hostremoved, INPUT_CHECK.out.fasta )
    }

    /*
        SUBWORKFLOW: PROFILING
    */

    PROFILING ( ch_reads_runmerged, DB_CHECK.out.dbs )
    ch_versions = ch_versions.mix( PROFILING.out.versions )

    /*
        SUBWORKFLOW: VISUALIZATION_KRONA
    */
    if ( params.run_krona ) {
        VISUALIZATION_KRONA ( PROFILING.out.classifications, PROFILING.out.profiles, DB_CHECK.out.dbs )
        ch_versions = ch_versions.mix( VISUALIZATION_KRONA.out.versions )
    }

    /*
        SUBWORKFLOW: PROFILING STANDARDISATION
    */
    if ( params.run_profile_standardisation ) {
        STANDARDISATION_PROFILES ( PROFILING.out.classifications, PROFILING.out.profiles, DB_CHECK.out.dbs, PROFILING.out.motus_version )
        ch_versions = ch_versions.mix( STANDARDISATION_PROFILES.out.versions )
    }

    /*
        MODULE: MultiQC
    */

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    workflow_summary    = WorkflowTaxprofiler.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowTaxprofiler.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    if ( params.preprocessing_qc_tool == 'falco' ) {
        // only mix in files actually used by MultiQC
        ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.txt
                            .map { meta, reports -> reports }
                            .flatten()
                            .filter { path -> path.name.endsWith('_data.txt')}
                            .ifEmpty([]))
    } else {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }


    if (params.perform_shortread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_longread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( LONGREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp'){
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_COMPLEXITYFILTERING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_shortread_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_HOSTREMOVAL.out.mqc.collect{it[1]}.ifEmpty([]))
    }

    if (params.perform_longread_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_HOSTREMOVAL.out.mqc.collect{it[1]}.ifEmpty([]))
    }

    ch_multiqc_files = ch_multiqc_files.mix( PROFILING.out.mqc.collect{it[1]}.ifEmpty([]) )

    if ( params.run_profile_standardisation ) {
        ch_multiqc_files = ch_multiqc_files.mix( STANDARDISATION_PROFILES.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
