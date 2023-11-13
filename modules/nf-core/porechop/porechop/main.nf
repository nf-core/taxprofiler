process PORECHOP_PORECHOP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::porechop=0.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h7cff6ad_2' :
        'biocontainers/porechop:0.2.4--py39h7cff6ad_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_porechopped.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## To ensure ID matches rest of pipeline based on meta.id rather than input file name
    
    [[ -f ${prefix}.fastq.gz   ]] || ln -s $reads ${prefix}.fastq.gz

    porechop \\
        -i ${prefix}.fastq.gz \\
        -t $task.cpus \\
        $args \\
        -o ${prefix}_porechopped.fastq.gz \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$( porechop --version )
    END_VERSIONS
    """
}
