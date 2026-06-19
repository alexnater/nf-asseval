process BLOBTOOLS_SNAIL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blobtoolkit:4.5.3--pyhdfd78af_0' :
        'biocontainers/blobtoolkit:4.5.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(full_table)

    output:
    tuple val(meta), path("*.snail.svg")      , optional: true, emit: snail
    tuple val("${task.process}"), val("blobtoolkit"), eval("blobtools --version | sed 's/.*v//'"), topic: versions, emit: versions_blobtoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix        = task.ext.prefix ?: "${meta.id}"
    def args      = task.ext.args ?: ''
    def args2     = task.ext.args2 ?: ''
    def busco_arg = full_table ? "--busco ${full_table}"  : ""

    """
    blobtools create \\
        --fasta ${fasta} \\
        $busco_arg \\
        $args \\
        ${prefix}

    blobtools view \\
        --plot \\
        --format svg \\
        --view snail \\
        $args2 \\
        ${prefix}
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.snail.svg
    """
}
