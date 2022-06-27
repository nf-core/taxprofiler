process EXTRACT_READLENGTH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::bash=5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv2/biocontainers_v1.2.0_cv2.img' :
        'biocontainers/biocontainers:v1.2.0_cv2' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), stdout, emit: read_length

    script:
    """
    tail -n 1 $tsv | cut -f 7 | tr -d '\n'
    """
}
