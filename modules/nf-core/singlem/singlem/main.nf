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

    echo "Metapackage root directory: \$MP_DIR"
    echo "Directory tree (max depth 3):"
    find "\$MP_DIR" -maxdepth 3 -print

    # Find directory that looks like a SingleM metapackage (*.smpkg.zb)
    MP_CAND=\$(find "\$MP_DIR" -maxdepth 3 -type d -name '*.smpkg.zb' | head -n 1)

    # Fallback: maybe MP_DIR itself is the metapackage
    if [ -z "\$MP_CAND" ]; then
        if [ -f "\$MP_DIR/CONTENTS.json" ] && [ -d "\$MP_DIR/payload_directory" ]; then
            MP_CAND="\$MP_DIR"
        fi
    fi

    echo "Using SingleM metapackage path: \$MP_CAND"

    if [ -z "\$MP_CAND" ]; then
        echo "ERROR: Could not locate SingleM metapackage directory under \$MP_DIR" >&2
        exit 1
    fi

    export SINGLEM_METAPACKAGE_PATH="\$MP_CAND"

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