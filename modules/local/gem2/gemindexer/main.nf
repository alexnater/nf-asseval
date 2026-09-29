process GEM2_GEMINDEXER {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/239a50bbee5154eb584f8d7c46d9c4f6d5a563ac8b3f9b3f883fe308cbc7ef58/data':
        'community.wave.seqera.io/library/gem2_ucsc-bedgraphtobigwig:d215a32896887bf0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gem"), emit: index
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('gem2'), val('20200110'), emit: versions_gem2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gem-indexer \\
        -i ${fasta} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        --mm-tmp-prefix ./tmp \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20200110' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.gem
    touch ${prefix}.log
    """
}
