process SINGLEM_PIPE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://wwood/singlem:0.20.3' :
        'docker://wwood/singlem:0.20.3' }"

    input:
    tuple val(meta), path(reads), path(metapackage_dir)

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
    MP_DIR="$metapackage_dir"
    MP_FILE=\$(find "\$MP_DIR" -maxdepth 2 -type f -name '*.smpkg.zb' | head -n 1)

    echo "Metapackage directory: \$MP_DIR"
    echo "Using metapackage file: \$MP_FILE"

    if [ -z "\$MP_FILE" ]; then
        echo "ERROR: No .smpkg.zb metapackage file found under \$MP_DIR" >&2
        exit 1
    fi

    export SINGLEM_METAPACKAGE_PATH="\$MP_FILE"

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