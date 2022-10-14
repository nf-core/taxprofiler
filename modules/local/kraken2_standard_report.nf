process KRAKEN2_STANDARD_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? 'conda-forge::sed=4.8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv2/biocontainers_v1.2.0_cv2.img' :
        'biocontainers/biocontainers:v1.2.0_cv2' }"

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

