//
// Subworkflow with functionality specific to the nf-core/pipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    databases         //  string: Path to databases

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --databases databases.csv --outdir <OUTDIR> --run_kraken2"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions

// Create channel from databases file provided through params.databases
    Channel
        .fromSamplesheet("databases")
        .set {ch_databases}

    emit:
    databases   = ch_databases
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs, fasta) = input[1..3]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs, fasta ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def text_seq_qc = [
        "Sequencing quality control was carried out with:",
        params.preprocessing_qc_tool == "falco" ? "Falco (de Sena Brandine and Smith 2021)." : "FastQC (Andrews 2010)."
    ].join(' ').trim()

    def text_shortread_qc = [
        "Short read preprocessing was performed with:",
        params.shortread_qc_tool == "adapterremoval" ? "AdapterRemoval (Schubert et al. 2016)." : "",
        params.shortread_qc_tool == "fastp" ? "fastp (Chen et al. 2018)." : "",
    ].join(' ').trim()

    def text_longread_qc = [
        "Long read preprocessing was performed with:",
        !params.longread_qc_skipadaptertrim ? "Porechop (Wick et al. 2017)," : "",
        !params.longread_qc_skipqualityfilter ? "Filtlong (Wick 2021)," : "",
        "."
    ].join(' ').trim()

    def text_shortreadcomplexity = [
        "Low-complexity sequence filtering was carried out with:",
        params.shortread_complexityfilter_tool == "bbduk" ? "BBDuk (Bushnell 2022)." : "",
        params.shortread_complexityfilter_tool == "prinseqplusplus" ? "PRINSEQ++ (Cantu et al. 2019)." : "",
        params.shortread_complexityfilter_tool == "fastp" ? "fastp (Chen et al. 2018)." : "",
    ].join(' ').trim()

    def text_shortreadhostremoval = [
        "Host read removal was performed for short reads with Bowtie2 (Langmead and Salzberg 2012) and SAMtools (Danecek et al. 2021)."
    ].join(' ').trim()

    def text_longreadhostremoval = [
        "Host read removal was performed for long reads with minimap2 (Li et al. 2018) and SAMtools (Danecek et al. 2021)."
    ].join(' ').trim()


    def text_classification = [
            "Taxonomic classification or profiling was carried out with:",
            params.run_bracken    ? "Bracken (Lu et al. 2017)," : "",
            params.run_kraken2    ? "Kraken2 (Wood et al. 2019)," : "",
            params.run_krakenuniq ? "KrakenUniq (Breitwieser et al. 2018)," : "",
            params.run_metaphlan  ? "MetaPhlAn (Blanco-Míguez et al. 2023)," : "",
            params.run_malt       ? "MALT (Vågene et al. 2018) and MEGAN6 CE (Huson et al. 2016)," : "",
            params.run_diamond    ? "DIAMOND (Buchfink et al. 2015)," : "",
            params.run_centrifuge ? "Centrifuge (Kim et al. 2016)," : "",
            params.run_kaiju      ? "Kaiju (Menzel et al. 2016)," : "",
            params.run_motus      ? "mOTUs (Ruscheweyh et al. 2022)," : "",
            params.run_ganon      ? "ganon (Piro et al. 2020)" : "",
            params.run_kmcp       ? "KMCP (Shen et al. 2023)" : "",
        "."
    ].join(' ').trim()

    def text_visualisation = [
            "Visualisation of results, where supported, was performed with Krona (Ondov et al. 2011)."
    ].join(' ').trim()

    def text_postprocessing = [
        "Standardisation of taxonomic profiles was carried out with TAXPASTA (Beber et al. 2023).",
    ].join(' ').trim()

    def citation_text = [
        "Tools used in the workflow included:",
        text_seq_qc,
        params.perform_shortread_qc               ? text_shortread_qc : "",
        params.perform_longread_qc                ? text_longread_qc : "",
        params.perform_shortread_complexityfilter ? text_shortreadcomplexity : "",
        params.perform_shortread_hostremoval      ? text_shortreadhostremoval : "",
        params.perform_longread_hostremoval       ? text_longreadhostremoval : "",
        text_classification,
        params.run_krona                          ? text_visualisation : "",
        params.run_profile_standardisation        ? text_postprocessing : "",
        "Pipeline results statistics were summarised with MultiQC (Ewels et al. 2016)."
    ].join(' ').trim().replaceAll("[,|.] +\\.", ".")

    return citation_text
}

def toolBibliographyText() {
    def text_seq_qc = [
        params.preprocessing_qc_tool == "falco"  ? "<li>de Sena Brandine, G., & Smith, A. D. (2021). Falco: high-speed FastQC emulation for quality control of sequencing data. F1000Research, 8(1874), 1874.  <a href=\"https://doi.org/10.12688/f1000research.21142.2\">10.12688/f1000research.21142.2</li>" : "",
        params.preprocessing_qc_tool == "fastqc" ? "<li>Andrews S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data, URL: <a href=\"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/\">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></li>" : "",
    ].join(' ').trim()


    def text_shortread_qc = [
        params.shortread_qc_tool == "adapterremoval" ? "<li>Schubert, M., Lindgreen, S., & Orlando, L. (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 9, 88. <a href=\"https://doi.org/10.1186/s13104-016-1900-2\">10.1186/s13104-016-1900-2</a></li>" : "",
    ].join(' ').trim()

    def text_longread_qc = [
        !params.longread_qc_skipadaptertrim   ? "<li>Wick, R. R., Judd, L. M., Gorrie, C. L., & Holt, K. E. (2017). Completing bacterial genome assemblies with multiplex MinION sequencing. Microbial Genomics, 3(10), e000132. <a href=\"https://doi.org/10.1099/mgen.0.000132\">10.1099/mgen.0.000132</a></li>" : "",
        !params.longread_qc_skipqualityfilter ? "<li>Wick R. (2021) Filtlong, URL:  <a href=\"https://github.com/rrwick/Filtlong\">https://github.com/rrwick/Filtlong</a></li>" : ""
    ].join(' ').trim()

    def text_shortreadcomplexity = [
        params.shortread_complexityfilter_tool == "bbduk" ? "<li>Bushnell B. (2022) BBMap, URL:  <a href=\"http://sourceforge.net/projects/bbmap/\">http://sourceforge.net/projects/bbmap/</a></li>" : "",
        params.shortread_complexityfilter_tool == "prinseqplusplus" ? "<li>Cantu, V. A., Sadural, J., & Edwards, R. (2019). PRINSEQ++, a multi-threaded tool for fast and efficient quality control and preprocessing of sequencing datasets (e27553v1). PeerJ Preprints. <a href=\"https://doi.org/10.7287/peerj.preprints.27553v1\">10.7287/peerj.preprints.27553v1</a></li>" : "",
    ].join(' ').trim()

    def text_shortreadhostremoval = [
        "<li>Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. <a href=\"https://doi.org/10.1038/nmeth.1923\">10.1038/nmeth.1923</a></li>",
    ].join(' ').trim()

    def text_longreadhostremoval = [
        "<li>Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics , 34(18), 3094–3100. <a href=\"https://doi.org/10.1093/bioinformatics/bty191\">10.1093/bioinformatics/bty191</a></li>",
    ].join(' ').trim()


    def text_classification = [
        params.run_bracken    ? "<li>Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ. Computer Science, 3(e104), e104. <a href=\"https://doi.org/10.7717/peerj-cs.104\">10.7717/peerj-cs.104</a></li>" : "",
        params.run_kraken2    ? "<li>Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1), 257.  <a href=\"https://doi.org/10.1186/s13059-019-1891-0\">10.1186/s13059-019-1891-0</a></li>" : "",
        params.run_krakenuniq ? "<li>Breitwieser, F. P., Baker, D. N., & Salzberg, S. L. (2018). KrakenUniq: confident and fast metagenomics classification using unique k-mer counts. Genome Biology, 19(1), 198.  <a href=\"https://doi.org/10.1186/s13059-018-1568-0\">10.1186/s13059-018-1568-0</a></li>" : "",
        params.run_metaphlan  ? "<li>Blanco-Míguez, A., Beghini, F., Cumbo, F., McIver, L. J., Thompson, K. N., Zolfo, M., Manghi, P., Dubois, L., Huang, K. D., Thomas, A. M., Nickols, W. A., Piccinno, G., Piperni, E., Punčochář, M., Valles-Colomer, M., Tett, A., Giordano, F., Davies, R., Wolf, J., … Segata, N. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. Nature Biotechnology, 1–12. <a href=\"https://doi.org/10.1038/s41587-023-01688-w\">10.1038/s41587-023-01688-w</a></li>" : "",
        params.run_malt       ? "<li>Vågene, Å. J., Herbig, A., Campana, M. G., Robles García, N. M., Warinner, C., Sabin, S., Spyrou, M. A., Andrades Valtueña, A., Huson, D., Tuross, N., Bos, K. I., & Krause, J. (2018). Salmonella enterica genomes from victims of a major sixteenth-century epidemic in Mexico. Nature Ecology & Evolution, 2(3), 520–528.  <a href=\"https://doi.org/10.1038/s41559-017-0446-6\">10.1038/s41559-017-0446-6</a></li>" : "",
        params.run_malt       ? "<li>Huson, D. H., Beier, S., Flade, I., Górska, A., El-Hadidi, M., Mitra, S., Ruscheweyh, H.-J., & Tappu, R. (2016). MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data. PLoS Computational Biology, 12(6), e1004957. <a href=\"https://doi.org/10.1371/journal.pcbi.1004957\">10.1371/journal.pcbi.1004957</a></li>" : "",
        params.run_diamond    ? "<li>Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59–60. <a href=\"https://doi.org/10.1038/nmeth.3176\">10.1038/nmeth.3176</a></li>" : "",
        params.run_centrifuge ? "<li>Kim, D., Song, L., Breitwieser, F. P., & Salzberg, S. L. (2016). Centrifuge: rapid and sensitive classification of metagenomic sequences. Genome Research, 26(12), 1721–1729.  <a href=\"https://doi.org/10.1101/gr.210641.116\">10.1101/gr.210641.116</a></li>" : "",
        params.run_kaiju      ? "<li>Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nature Communications, 7, 11257. <a href=\"https://doi.org/10.1038/ncomms11257\">10.1038/ncomms11257</a></li>" : "",
        params.run_motus      ? "<li>Ruscheweyh, H.-J., Milanese, A., Paoli, L., Karcher, N., Clayssen, Q., Keller, M. I., Wirbel, J., Bork, P., Mende, D. R., Zeller, G., & Sunagawa, S. (2022). Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments. Microbiome, 10(1), 212. <a href=\"https://doi.org/10.1186/s40168-022-01410-z\">10.1186/s40168-022-01410-z</a></li>" : "",
        params.run_ganon      ? "<li>Piro, V. C., Dadi, T. H., Seiler, E., Reinert, K., & Renard, B. Y. (2020). Ganon: Precise metagenomics classification against large and up-to-date sets of reference sequences. Bioinformatics (Oxford, England), 36(Suppl_1), i12–i20. <a href=\"https://doi.org/10.1093/bioinformatics/btaa458\">10.1093/bioinformatics/btaa458</a></li>" : "",
        params.run_kmcp       ? "<li>Shen, W., Xiang, H., Huang, T., Tang, H., Peng, M., Cai, D., Hu, P., & Ren, H. (2023). KMCP: accurate metagenomic profiling of both prokaryotic and viral populations by pseudo-mapping. Bioinformatics (Oxford, England), 39(1). <a href=\"https://doi.org/10.1093/bioinformatics/btac845\">10.1093/bioinformatics/btac845</a></li>" : "",
    ].join(' ').trim()

    def text_visualisation = [
        "<li>Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12(1), 385.  <a href=\"https://doi.org/10.1186/1471-2105-12-385\">10.1186/1471-2105-12-385</a></li>"
    ].join(' ').trim()

    def text_postprocessing = [
        "<li>Beber, M. E., Borry, M., Stamouli, S., & Fellows Yates, J. A. (2023). TAXPASTA: TAXonomic Profile Aggregation and STAndardisation. Journal of Open Source Software, 8(87), 5627.  <a href=\"https://doi.org/10.21105/joss.05627\">10.21105/joss.05627</a></li>",
    ].join(' ').trim()

    def text_extras = [
        // fastp shortread qc / complexity filtering
        ( params.perform_shortread_qc && params.shortread_qc_tool == "fastp" ) || ( params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool == "fastp" ) ? "<li>Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics , 34(17), i884–i890. <a href=\"https://doi.org/10.1093/bioinformatics/bty560\">10.1093/bioinformatics/bty560</a></li>" : "",
        // samtools long / short hostremoval
        params.perform_shortread_hostremoval || params.perform_longread_hostremoval ? "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2). <a href=\"https://doi.org/10.1093/gigascience/giab008\">10.1093/gigascience/giab008</a></li>" : "",
    ].join(' ').trim()

    def reference_text = [
        text_seq_qc,
        params.perform_shortread_qc               ? text_shortread_qc : "",
        params.perform_longread_qc                ? text_longread_qc : "",
        params.perform_shortread_complexityfilter ? text_shortreadcomplexity : "",
        params.perform_shortread_hostremoval      ? text_shortreadhostremoval : "",
        params.perform_longread_hostremoval       ? text_longreadhostremoval : "",
        text_extras,
        text_classification,
        params.run_krona                          ? text_visualisation : "",
        params.run_profile_standardisation        ? text_postprocessing : "",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. <a href=\"https:/doi.org/10.1093/bioinformatics/btw354\">10.1093/bioinformatics/btw354.</a></li>"
    ].join(' ').trim().replaceAll("[,|.] +\\.", ".")

    return reference_text
}

def methodsDescriptionText( mqc_methods_yaml ) {

    def meta = [:]

    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>Stamouli et al. 2023</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
