process SMUDGEPLOT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45ca1be680749bc69ce0557205ffc6e5778e5032d3134ff8966258de51aa46cc/data' :
        'community.wave.seqera.io/library/fastk_smudgeplot:c597618c48a0dba3' }"

    input:
    tuple val(meta), path(ktab), path(data)
    val(cutoff)

    output:
    tuple val(meta), path("*.smu")            , emit: smu
    tuple val(meta), path("*.sma")            , emit: sma
    tuple val(meta), path("*.tsv")            , emit: report
    tuple val(meta), path("*.png")            , emit: plots
    tuple val(meta), path("*.txt")            , emit: txt
    tuple val("${task.process}"), val('smudgeplot'), eval('smudgeplot.py --version 2>&1 | sed "s/^.*smudgeplot //"'), emit: versions_smudgeplot, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p tmp
    export MPLCONFIGDIR=./tmp

    smudgeplot \\
        hetmers \\
        -t $task.cpus \\
        -tmp tmp \\
        -L $cutoff \\
        -o $prefix \\
        $args \\
        $ktab

    smudgeplot \\
        all \\
        -o $prefix \\
        $args2 \\
        ${prefix}.smu

    rm -r tmp
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.smu
    touch ${prefix}.smudge_report.tsv
    touch ${prefix}_centralities.txt
    touch ${prefix}_centralities.png
    touch ${prefix}_smudgeplot.png
    touch ${prefix}_smudgeplot_log10.png
    """
}
