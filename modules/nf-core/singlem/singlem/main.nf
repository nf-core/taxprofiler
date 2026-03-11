process SINGLEM_PIPE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://wwood/singlem:0.20.3' :
        'docker://wwood/singlem:0.20.3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.profile.tsv'), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_args = meta.single_end ? "-1 ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    singlem data
    export SINGLEM_METAPACKAGE_PATH=~/.singlem
    singlem pipe \\
        $input_args \\
        -p ${prefix}.profile.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1)
    END_VERSIONS
    """
}