//
// This file holds several functions specific to the workflow/taxprofiler.nf in the nf-core/taxprofiler pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowTaxprofiler {

    //
    // Check and validate parameters
    //

    public static void initialise(params, log) {

        genomeExistsError(params, log)
        //if (!params.fasta) {
        //    Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        //}
    }


    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    ///
    /// Automatic publication methods text generation
    ///

    public static String toolCitationText(params) {

        def citation_text = [
                "Tools used in the workflow included",
                params["preprocessing_qc_tool"] == "falco" ? "Falco (de Sena Brandine and Smith 2021)." : "FastQC (Andrews 2010).", // TODO OR FALCO

                params["perform_shortread_qc"] ? ". Short read preprocessing was carried out with" : "",
                params["perform_shortread_qc"] && params["shortread_qc_tool"] == "adapterremoval" ? "AdapterRemoval (Schubert et al. 2016)." : "",
                params["perform_shortread_qc"] && params["shortread_qc_tool"] == "fastp" ? "fastp (Chen et al. 2018)." : "",

                params["perform_longread_qc"] ? ". Long read preprocessing was carried out with" : "",
                params["perform_longread_qc"] && !params["longread_qc_skipadaptertrim"] ? "Porechop (Wick 2018)," : "",
                params["perform_longread_qc"] && !params["longread_qc_skipqualityfilter"] ? "Filtlong (Wick 2021)," : "",

                params["perform_shortread_complexityfilter"] ? ". Complexity filtering was performed using" : "",
                params["perform_shortread_complexityfilter"] && params["shortread_complexityfilter_tool"] == "bbduk" ? "BBDuk (Bushnell 2022)," : "",
                params["perform_shortread_complexityfilter"] && params["shortread_complexityfilter_tool"] == "prinseqplusplus" ? "PRINSEQ++ (Cantu et al. 2019)," : "",
                params["perform_shortread_complexityfilter"] && params["shortread_complexityfilter_tool"] == "fastp" ? "fastp (Chen et al. 2018)," : "",

                params["perform_shortread_hostremoval"] ? ". Host read removal was carried out for short reads with Bowtie2 (Langmead and Salzberg 2012)" : "",
                params["perform_longread_hostremoval"] ? ". Host read removal was carried out for long reads with minimap2 (Li et al. 2018)" : "",
                params["perform_shortread_hostremoval"] || params["perform_longread_hostremoval"] ? "and SAMtools (Danecek et al. 2021)." : "",

                ". Taxonomic classification or profiling was performed with",
                params["run_bracken"] ? "Bracken (Lu et al. 2017)," : "",
                params["run_kraken2"] ? "Kraken2 (Wood et al. 2019)," : "",
                params["run_krakenuniq"] ? "KrakenUniq (Breitwieser et al. 2018)," : "",
                params["run_metaphlan3"] ? "MetaPhlAn3 (Beghini et al. 2021)," : "",
                params["run_malt"] ? "MALT (Vågene et al. 2018) and MEGAN6 CE (Huson et al. 2016)," : "",
                params["run_diamond"] ? "DIAMOND (Buchfink et al. 2015)," : "",
                params["run_centrifuge"] ? "Centrifuge (Kim et al. 2016)," : "",
                params["run_kaiju"] ? "Kaiju (Menzel et al. 2016)," : "",
                params["run_motus"] ? "mOTUs (Ruscheweyh et al. 2022)," : "",

                params["run_krona"] ? ". Results visualisation for some tools were displayed with Krona (Ondov et al. 2011)." : "",

                ". Pipeline results statistics were summarised with MultiQC (Ewels et al. 2016)."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // TODO consider how to do the same for the references themselves, include in the same if/else statements somehow?
        def reference_text = [
                params["preprocessing_qc_tool"] == "falco" ? "<li>de Sena Brandine G and Smith AD. Falco: high-speed FastQC emulation for quality control of sequencing data. F1000Research 2021, 8:1874</li>" : "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>", // TODO OR FALCO

                params["perform_shortread_qc"] && params["shortread_qc_tool"] == "adapterremoval" ? "<li>Schubert, Mikkel, Stinus Lindgreen, and Ludovic Orlando. 2016. AdapterRemoval v2: Rapid Adapter Trimming, Identification, and Read Merging. BMC Research Notes 9 (February): 88. doi:10.1186/s13104-016-1900-2.</li>" : "",
                params["perform_shortread_qc"] && params["shortread_qc_tool"] == "fastp" ? "<li>Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. 2018. Fastp: An Ultra-Fast All-in-One FASTQ Preprocessor. Bioinformatics 34 (17): i884-90. 10.1093/bioinformatics/bty560.</li>" : "",

                params["perform_longread_qc"] && !params["longread_qc_skipadaptertrim"] ? "<li>Wick R (2018) Porechop, URL: https://github.com/rrwick/Porechop</li>" : "",
                params["perform_longread_qc"] && !params["longread_qc_skipqualityfilter"] ? "<li>Wick R (2021) Filtlong, URL: https://github.com/rrwick/Filtlong</li>" : "",

                params["perform_shortread_complexityfilter"] && params["shortread_complexityfilter_tool"] == "bbduk" ? "<li>Bushnell B (2022) BBMap, URL: http://sourceforge.net/projects/bbmap/</li>" : "",
                params["perform_shortread_complexityfilter"] && params["shortread_complexityfilter_tool"] == "prinseqplusplus" ? "<li>Cantu, Vito Adrian, Jeffrey Sadural, and Robert Edwards. 2019. PRINSEQ++, a Multi-Threaded Tool for Fast and Efficient Quality Control and Preprocessing of Sequencing Datasets. e27553v1. PeerJ Preprints. doi: 10.7287/peerj.preprints.27553v1.</li>" : "",
                params["perform_shortread_complexityfilter"] && params["shortread_complexityfilter_tool"] == "fastp" ? "<li>Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. 2018. Fastp: An Ultra-Fast All-in-One FASTQ Preprocessor. Bioinformatics 34 (17): i884-90. 10.1093/bioinformatics/bty560.</li>" : "",

                params["perform_shortread_hostremoval"] ? "<li>Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. doi: 10.1038/nmeth.1923</li>" : "",
                params["perform_longread_hostremoval"] ? "<li>Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics , 34(18), 3094–3100. doi: 10.1093/bioinformatics/bty191</li>" : "",
                params["perform_shortread_hostremoval"] || params["perform_longread_hostremoval"] ? "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2). doi: 10.1093/gigascience/giab008</li>" : "",

                params["run_bracken"] ? "<li>Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: Estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. doi: 10.7717/peerj-cs.104</li>" : "",
                params["run_kraken2"] ? "<li>Wood, Derrick E., Jennifer Lu, and Ben Langmead. 2019. Improved Metagenomic Analysis with Kraken 2. Genome Biology 20 (1): 257. doi: 10.1186/s13059-019-1891-0.</li>" : "",
                params["run_krakenuniq"] ? "<li>Breitwieser, Florian P., Daniel N. Baker, and Steven L. Salzberg. 2018. KrakenUniq: confident and fast metagenomics classification using unique k-mer counts. Genome Biology 19 (1): 198. doi: 10.1186/s13059-018-1568-0</li>" : "",
                params["run_metaphlan3"] ? "<li>Beghini, Francesco, Lauren J McIver, Aitor Blanco-Míguez, Leonard Dubois, Francesco Asnicar, Sagun Maharjan, Ana Mailyan, et al. 2021. “Integrating Taxonomic, Functional, and Strain-Level Profiling of Diverse Microbial Communities with BioBakery 3.” Edited by Peter Turnbaugh, Eduardo Franco, and C Titus Brown. ELife 10 (May): e65088. doi: 10.7554/eLife.65088</li>" : "",
                params["run_malt"] ? "<li>Vågene, Åshild J., Alexander Herbig, Michael G. Campana, Nelly M. Robles García, Christina Warinner, Susanna Sabin, Maria A. Spyrou, et al. 2018. Salmonella Enterica Genomes from Victims of a Major Sixteenth-Century Epidemic in Mexico. Nature Ecology & Evolution 2 (3): 520-28. doi: 10.1038/s41559-017-0446-6.</li>" : "",
                params["run_malt"] ? "<li>Huson, Daniel H., Sina Beier, Isabell Flade, Anna Górska, Mohamed El-Hadidi, Suparna Mitra, Hans-Joachim Ruscheweyh, and Rewati Tappu. 2016. “MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data.” PLoS Computational Biology 12 (6): e1004957. doi: 10.1371/journal.pcbi.1004957.</li>" : "",
                params["run_diamond"] ? "<li>Buchfink, Benjamin, Chao Xie, and Daniel H. Huson. 2015. “Fast and Sensitive Protein Alignment Using DIAMOND.” Nature Methods 12 (1): 59-60. doi: 10.1038/nmeth.3176.</li>" : "",
                params["run_centrifuge"] ? "<li>Kim, Daehwan, Li Song, Florian P. Breitwieser, and Steven L. Salzberg. 2016. “Centrifuge: Rapid and Sensitive Classification of Metagenomic Sequences.” Genome Research 26 (12): 1721-29. doi: 10.1101/gr.210641.116.</li>" : "",
                params["run_kaiju"] ? "<li>Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nature Communications, 7, 11257. doi: 10.1038/ncomms11257</li>" : "",
                params["run_motus"] ? "<li>Ruscheweyh, H.-J., Milanese, A., Paoli, L., Karcher, N., Clayssen, Q., Keller, M. I., Wirbel, J., Bork, P., Mende, D. R., Zeller, G., & Sunagawa, S. (2022). Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments. Microbiome, 10(1), 212. doi: 10.1186/s40168-022-01410-z</li>" : "",
                params["run_krona"] ? "<li>Ondov, Brian D., Nicholas H. Bergman, and Adam M. Phillippy. 2011. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12 (1): 385. doi: 10.1186/1471-2105-12-385.</li>" : "",

                "<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]

        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""
        /*
        TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        */
        meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        meta["tool_bibliography"] = toolBibliographyText(params)

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }
}
