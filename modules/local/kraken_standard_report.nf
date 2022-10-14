process KRAKEN_STANDARD_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'conda-forge::sed=4.8' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv2/biocontainers_v1.2.0_cv2.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv2'
    }

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path(result), emit: report

    script:
    result = "${report.baseName}_standardized.kraken2.report.txt"
    """
    cut -f1-3,6-8 "${report}" > "${result}"
    """
}

