process NONPAREIL_TOJSON {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.4.1--r42h4ac6f70_4':
        'biocontainers/nonpareil:3.4.1--r42h4ac6f70_4' }"

    input:
    tuple val(meta), path(npo)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    ## Custom code from Filipe G. Vieira
    ## As described on the MultiQC website
    ## https://multiqc.info/modules/nonpareil/

    base::library("jsonlite")
    base::message("Exporting model as JSON")

    export_curve <- function(object){
    # Extract variables
    n <- names(attributes(object))[c(1:12,21:29)]
    x <- sapply(n, function(v) attr(object,v))
    names(x) <- n
    # Extract vectors
    n <- names(attributes(object))[13:20]
    y <- lapply(n, function(v) attr(object,v))
    names(y) <- n
    curve_json <- c(x, y)

    # Add model
    if (object\$has.model) {
        # https://github.com/lmrodriguezr/nonpareil/blob/162f1697ab1a21128e1857dd87fa93011e30c1ba/utils/Nonpareil/R/Nonpareil.R#L330-L332
        x_min <- 1e3
        x_max <- signif(tail(attr(object,"x.adj"), n=1)*10, 1)
        x.model <- exp(seq(log(x_min), log(x_max), length.out=1e3))
        y.model <- predict(object, lr=x.model)
        curve_json <- append(curve_json, list(x.model=x.model))
        curve_json <- append(curve_json, list(y.model=y.model))
    }

    base::print(curve_json)
    curve_json
    }

    export_set <- function(object){
    y <- lapply(object\$np.curves, "export_curve")
    names(y) <- sapply(object\$np.curves, function(n) n\$label)
    jsonlite::prettify(toJSON(y, auto_unbox=TRUE))
    }

    curves <- Nonpareil.curve($npo)

    y <- export_set(curves)
    write(y, "${prefix}.json")

    # version export
    f <- file("versions.yml","w")
    nonpareil_version = sessionInfo()\$otherPkgs\$Nonpareil\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    nonpareil:", nonpareil_version), f)
    close(f)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    file.create("${process}.json")

    # version export
    f <- file("versions.yml","w")
    nonpareil_version = sessionInfo()\$otherPkgs\$Nonpareil\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    nonpareil:", nonpareil_version), f)
    close(f)
    """
}
