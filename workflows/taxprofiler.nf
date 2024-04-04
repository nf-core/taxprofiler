/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_taxprofiler_pipeline'
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

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

if ( [params.taxpasta_add_name, params.taxpasta_add_rank, params.taxpasta_add_lineage, params.taxpasta_add_lineage, params.taxpasta_add_idlineage, params.taxpasta_add_ranklineage].any() && !params.taxpasta_taxonomy_dir ) error('ERROR: [nf-core/taxprofiler] All --taxpasta_add_* parameters require a taxonomy supplied to --taxpasta_taxonomy_dir. However the latter parameter was not detected. Please check input.')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

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
include { UNTAR                       } from '../modules/nf-core/untar/main'
include { FALCO                       } from '../modules/nf-core/falco/main'
include { CAT_FASTQ as MERGE_RUNS     } from '../modules/nf-core/cat/fastq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TAXPROFILER {

    adapterlist = params.shortread_qc_adapterlist ? file(params.shortread_qc_adapterlist) : []

    if ( params.shortread_qc_adapterlist ) {
        if ( params.shortread_qc_tool == 'adapterremoval' && !(adapterlist.extension == 'txt') ) error "[nf-core/taxprofiler] ERROR: AdapterRemoval2 adapter list requires a `.txt` format and extension. Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
        if ( params.shortread_qc_tool == 'fastp' && !adapterlist.extension.matches(".*(fa|fasta|fna|fas)") ) error "[nf-core/taxprofiler] ERROR: fastp adapter list requires a `.fasta` format and extension (or fa, fas, fna). Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
    }

    take:
    samplesheet // channel: samplesheet read in from --input
    databases // channel: databases from --databases

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Validate input files and create separate channels for FASTQ, FASTA, and Nanopore data
    ch_input = samplesheet
        .map { meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            meta.run_accession = run_accession
            meta.instrument_platform = instrument_platform

            // Define single_end based on the conditions
            meta.single_end = ( fastq_1 && !fastq_2 && instrument_platform != 'OXFORD_NANOPORE' )

            // Define is_fasta based on the presence of fasta
            meta.is_fasta = fasta ? true : false

            if ( !meta.is_fasta && !fastq_1 ) {
                error("ERROR: Please check input samplesheet: entry `fastq_1` doesn't exist!")
            }
            if ( meta.instrument_platform == 'OXFORD_NANOPORE' && fastq_2 ) {
                error("Error: Please check input samplesheet: for Oxford Nanopore reads entry `fastq_2` should be empty!")
            }
            if ( meta.single_end && fastq_2 ) {
                error("Error: Please check input samplesheet: for single-end reads entry `fastq_2` should be empty")
            }
            return [ meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ]
        }
        .branch { meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            fastq: meta.single_end || fastq_2
                return [ meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
            nanopore: instrument_platform == 'OXFORD_NANOPORE'
                meta.single_end = true
                return [ meta, [ fastq_1 ] ]
            fasta: meta.is_fasta
                meta.single_end = true
                return [ meta, [ fasta ] ]
        }

    // Merge ch_input.fastq and ch_input.nanopore into a single channel
    def ch_input_for_fastqc = ch_input.fastq.mix( ch_input.nanopore )

    // Validate and decompress databases
    ch_dbs_for_untar = databases
        .branch { db_meta, db_path ->
            untar: db_path.name.endsWith( ".tar.gz" )
            skip: true
        }
    // Filter the channel to untar only those databases for tools that are selected to be run by the user.
    // Also, to ensure only untar once per file, group together all databases of one file
    ch_inputdb_untar = ch_dbs_for_untar.untar
        .filter { db_meta, db_path ->
            params[ "run_${db_meta.tool}" ]
        }
        .groupTuple(by: 1)
        .map {
            meta, dbfile ->
                def new_meta = [ 'id': dbfile.baseName ] + [ 'meta': meta ]
            [new_meta , dbfile ]
        }

    // Untar the databases
    UNTAR ( ch_inputdb_untar )
    ch_versions = ch_versions.mix( UNTAR.out.versions.first() )

    // Spread out the untarred and shared databases
    ch_outputdb_from_untar = UNTAR.out.untar
        .map {
            meta, db ->
            [meta.meta, db]
        }
        .transpose(by: 0)

    ch_final_dbs = ch_dbs_for_untar.skip
                    .mix( ch_outputdb_from_untar  )
                    .map { db_meta, db ->
                        def corrected_db_params = db_meta.db_params ? [ db_params: db_meta.db_params ] : [ db_params: '' ]
                        [ db_meta + corrected_db_params, db ]
                    }

    /*
        MODULE: Run FastQC
    */


    if ( !params.skip_preprocessing_qc ) {
        if ( params.preprocessing_qc_tool == 'falco' ) {
            FALCO ( ch_input_for_fastqc )
            ch_versions = ch_versions.mix(FALCO.out.versions.first())
        } else {
            FASTQC ( ch_input_for_fastqc )
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        }
    }

    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */

    if ( params.perform_shortread_qc ) {
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING ( ch_input.fastq, adapterlist ).reads
        ch_versions = ch_versions.mix( SHORTREAD_PREPROCESSING.out.versions )
    } else {
        ch_shortreads_preprocessed = ch_input.fastq
    }

    if ( params.perform_longread_qc ) {
        ch_longreads_preprocessed = LONGREAD_PREPROCESSING ( ch_input.nanopore ).reads
                                        .map { it -> [ it[0], [it[1]] ] }
        ch_versions = ch_versions.mix( LONGREAD_PREPROCESSING.out.versions )
    } else {
        ch_longreads_preprocessed = ch_input.nanopore
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
                    def meta_new = meta - meta.subMap('run_accession')
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
            .mix( ch_input.fasta )

        ch_versions = ch_versions.mix(MERGE_RUNS.out.versions)

    } else {
        ch_reads_runmerged = ch_shortreads_hostremoved
            .mix( ch_longreads_hostremoved, ch_input.fasta )
    }

    /*
        SUBWORKFLOW: PROFILING
    */

    PROFILING ( ch_reads_runmerged, ch_final_dbs )
    ch_versions = ch_versions.mix( PROFILING.out.versions )

    /*
        SUBWORKFLOW: VISUALIZATION_KRONA
    */
    if ( params.run_krona ) {
        VISUALIZATION_KRONA ( PROFILING.out.classifications, PROFILING.out.profiles, ch_final_dbs )
        ch_versions = ch_versions.mix( VISUALIZATION_KRONA.out.versions )
    }

    /*
        SUBWORKFLOW: PROFILING STANDARDISATION
    */
    if ( params.run_profile_standardisation ) {
        STANDARDISATION_PROFILES ( PROFILING.out.classifications, PROFILING.out.profiles, ch_final_dbs, PROFILING.out.motus_version )
        ch_versions = ch_versions.mix( STANDARDISATION_PROFILES.out.versions )
    }

    /*
        MODULE: MultiQC
    */

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    if ( !params.skip_preprocessing_qc ) {
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

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
