process SINGLEM_PIPE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Use container if desired, optional
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/singlem:0.20.3--py39h6a87c0e_0':
        'biocontainers/singlem:0.20.3--py39h6a87c0e_0' }"

    input:
    tuple val(meta), path(reads)
    // path(db)  // optional if you want to provide a SingleM DB

    output:
    tuple val(meta), path('*.profile.tsv'), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // SingleM command arguments
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Determine input: paired vs single
    def input_args = meta.single_end ? "-1 ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    singlem pipe \\
        $input_args \\
        -p ${prefix}.profile.tsv \\
        $args

    # Save software version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1)
    END_VERSIONS
    """
}