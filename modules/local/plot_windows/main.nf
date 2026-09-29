process PLOT_WINDOWS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/11/115554413e53171dfff5f5752a67708e1e4e837c2fd87e333e310076110e54b7/data' :
        'community.wave.seqera.io/library/r-ggplot2_r-optparse:6570bf84e2757fb9' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.winstats.pdf"), emit: pdf
    tuple val("${task.process}"), val('rscript'), eval("Rscript --version | sed 's/^.*version \\(.*\\) (.*/\\1/'"), emit: versions_rscript, topic: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_windowstats.r \\
        --bed $bed \\
        --outfile ${prefix}.winstats.pdf \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.winstats.pdf
    """
}