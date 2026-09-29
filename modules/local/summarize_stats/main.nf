process SUMMARIZE_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d459c2a6af06632b91dc8e73d86f6b7bad18913e0482c40452767df3e3375004/data' :
        'community.wave.seqera.io/library/pysam_numpy:beea5885b4c2cb2a' }"

    input:
    tuple val(meta), path(flagstats), path(depths), path(coverages), val(before), val(after)

    output:
    tuple val(meta), path("*.report.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    summarize_stats.py \\
        -f $flagstats \\
        -s $depths \\
        -c $coverages \\
        -b ${before.join(' ')} \\
        -a ${after.join(' ')} \\
        $args \\
        -o ${prefix}.report.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.report.tsv
    """
}