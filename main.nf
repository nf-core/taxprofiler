#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/taxprofiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/taxprofiler
    Website: https://nf-co.re/taxprofiler
    Slack  : https://nfcore.slack.com/channels/taxprofiler
----------------------------------------------------------------------------------------
*/

nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TAXPROFILER             } from './workflows/taxprofiler'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_taxprofiler_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_taxprofiler_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.databases,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_TAXPROFILER(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.databases,
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_TAXPROFILER.out.multiqc_report,
    )

    processed_reads = NFCORE_TAXPROFILER.out.processed_reads
        .map { meta, file ->
            def files = [file].flatten()
            [
                id: meta.id,
                run_accession: meta.run_accession,
                instrument_platform: meta.instrument_platform,
                fastq_1: !meta.is_fasta ? files[0] : "",
                fastq_2: !meta.single_end && !meta.is_fasta ? files[1] : "",
                fasta: meta.is_fasta ? files[0] : "",
            ]
        }
    classifications = NFCORE_TAXPROFILER.out.classifications
        .map { meta, file ->
            [
                id: meta.id,
                instrument_platform: meta.instrument_platform,
                tool: meta.tool,
                db_name: meta.db_name,
                file: file
            ]
        }
    profiles = NFCORE_TAXPROFILER.out.profiles
        .map { meta, file ->
            [
                id: meta.id,
                instrument_platform: meta.instrument_platform,
                tool: meta.tool,
                db_name: meta.db_name,
                file: file
            ]
        }
    profile_saved_reads = NFCORE_TAXPROFILER.out.profile_saved_reads
        .map { meta, file ->
            [
                id: meta.id,
                instrument_platform: meta.instrument_platform,
                tool: meta.tool,
                db_name: meta.db_name,
                file: file
            ]
        }

    taxpasta = Channel.empty()
    combined_report = Channel.empty()
    if ( params.run_profile_standardisation ) {
        taxpasta = taxpasta.mix(NFCORE_TAXPROFILER.out.taxpasta)
            .map { meta, file ->
                [
                    tool: meta.tool,
                    db_name: meta.db_name,
                    file: file
                ]
            }
        combined_report = combined_report.mix(NFCORE_TAXPROFILER.out.combined_report)
    }


    publish:
    processed_reads        = processed_reads
    classifications        = classifications
    profiles               = profiles
    bracken_kraken2_report = NFCORE_TAXPROFILER.out.bracken_kraken2_report
    centrifuge_report      = NFCORE_TAXPROFILER.out.centrifuge_report
    diamond_log            = NFCORE_TAXPROFILER.out.diamond_log
    ganon_outputs          = NFCORE_TAXPROFILER.out.ganon_outputs
    kmcp_search            = NFCORE_TAXPROFILER.out.kmcp_search
    metaphlan_outputs      = NFCORE_TAXPROFILER.out.metaphlan_outputs
    profile_saved_reads    = profile_saved_reads
    taxpasta               = taxpasta
    combined_report        = combined_report
}

output {
    'processed_reads' {
        path {rec ->
            if (rec.instrument_platform == "OXFORD_NANOPORE" && params.perform_longread_qc) {
                "${params.longread_filter_tool}"
            } else if (params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp') {
                "bbduk"
            } else {
                null
            }
        }
        index {
            path 'processed_reads.csv'
            header true
            sep ','
        }
    }
    'classifications' {
        path { rec ->
            if (rec.tool == "bracken") {
                "kraken2/${rec.db_name}_bracken"
            } else {
                "${rec.tool}/${rec.db_name}"
            }
        }
        index {
            path 'classifications.csv'
            header true
            sep ','
        }
    }
    'profiles' {
        path { rec ->
            if (rec.tool == "kraken2-bracken") {
                "kraken2/${rec.db_name}_bracken"
            } else {
                "${rec.tool}/${rec.db_name}"
            }
        }
        index {
            path 'profiles.csv'
            header true
            sep ','
        }
    }
    'bracken_kraken2_report' {
        path { meta, file ->
            "bracken/${meta.db_name}"
        }
    }
    'centrifuge_report' {
        path { meta, file ->
            "centrifuge/${meta.db_name}"
        }
    }
    'diamond_log' {
        path { meta, file ->
            "diamond/${meta.db_name}"

        }
    }
    'ganon_outputs' {
        path { meta, file ->
            "ganon/${meta.db_name}"
        }
    }
    'kmcp_search' {
        path { meta, file ->
            "kmcp/${meta.db_name}"
        }
    }
    'metaphlan_outputs' {
        path { meta, file ->
            "metaphlan/${meta.db_name}"

        }
    }
    'profile_saved_reads' {
        path { rec ->
            if (rec.tool == "bracken") {
                "kraken2/${rec.db_name}_bracken"
            } else {
                "${rec.tool}/${rec.db_name}"
            }
        }
    }
    'taxpasta' {
        path "taxpasta"
        enabled params.run_profile_standardisation
        index {
            path 'taxpasta.csv'
            header true
            sep ','
        }
    }
    'combined_report' {
        path { meta, file ->
            if (meta.tool == "kraken2-bracken") {
                "kraken2/"
            } else {
                "${meta.tool}"
            }
        }
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_TAXPROFILER {
    take:
    samplesheet // channel: samplesheet read in from --input
    databases   // channel: databases in from --databases

    main:

    //
    // WORKFLOW: Run pipeline
    //
    TAXPROFILER(
        samplesheet,
        databases,
    )

    emit:
    multiqc_report         = TAXPROFILER.out.multiqc_report // channel: /path/to/multiqc_report.html
    multiqc_data           = TAXPROFILER.out.multiqc_data
    multiqc_plots          = TAXPROFILER.out.multiqc_plots
    processed_reads        = TAXPROFILER.out.processed_reads
    classifications        = TAXPROFILER.out.classifications
    profiles               = TAXPROFILER.out.profiles
    bracken_kraken2_report = TAXPROFILER.out.bracken_kraken2_report
    centrifuge_report      = TAXPROFILER.out.centrifuge_report
    diamond_log            = TAXPROFILER.out.diamond_log
    ganon_outputs          = TAXPROFILER.out.ganon_outputs
    kmcp_search            = TAXPROFILER.out.kmcp_search
    metaphlan_outputs      = TAXPROFILER.out.metaphlan_outputs
    profile_saved_reads    = TAXPROFILER.out.profile_saved_reads
    taxpasta               = TAXPROFILER.out.taxpasta
    combined_report        = TAXPROFILER.out.combined_report
}
