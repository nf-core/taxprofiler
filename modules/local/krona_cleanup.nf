process KRONA_CLEANUP {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(krona, stageAs: 'uncleaned.krona.txt')

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Copy the file to a new name
    cp ${krona} ${prefix}.txt

    # Remove ugly 'x__' prefixes for each of the taxonomic levels
    LEVELS=(d k p c o f g s)
    for L in "\${LEVELS[@]}"; do
        sed -i "s/\${L}__//g" ${prefix}.txt
    done

    # Remove underscores that are standing in place of spaces
    sed -i "s/_/ /g" ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
