process KRAKEN2_STANDARD_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path(result), emit: report
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    result = "${prefix}_standardized.kraken2.report.txt"
    """
    cut -f1-3,6-8 '${report}' > '${result}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(echo \$(cut --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}

