process ENSURE_FASTQ_EXTENSION {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::bash=5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv2/biocontainers_v1.2.0_cv2.img' :
        'biocontainers/biocontainers:v1.2.0_cv2' }"


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads

    script:
    if (meta.single_end) {
        fastq = "${reads.baseName}.fastq.gz"
        """
        ln -s '${reads}' '${fastq}'
        """
    } else {
        first = "${reads[0].baseName}.fastq.gz"
        second = "${reads[1].baseName}.fastq.gz"
        """
        ln -s '${reads[0]}' '${first}'
        ln -s '${reads[1]}' '${second}'
        """
    }
}
