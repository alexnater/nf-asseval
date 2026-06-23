process GENERATE_EAR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b1dd9a5353bf5b3d502c6ea059b88c0378d57985142e45d40304897f8f43b8e3/data'
        : 'community.wave.seqera.io/library/python_pip_pyyaml_argparse:d26cead1d8e9b0cb'}"

    input:
    tuple val(meta) , path(template)
    tuple val(meta2), path(genomescope, stageAs: "data/kmer/*"), path(smudgeplot, stageAs: "data/kmer/*")
    tuple val(meta3), path('data/contig/merqury/*'), path('data/contig/stats.hap1.txt'), path('data/contig/stats.hap2.txt'), path('data/contig/busco.hap1.txt'), path('data/contig/busco.hap2.txt'), path('data/contig/snail.hap1.svg'), path('data/contig/snail.hap2.svg'), path('data/contig/blob.hap1.svg'), path('data/contig/blob.hap2.svg'), path('data/contig/pretext.hap1.png'), path('data/contig/pretext.hap2.png')
    tuple val(meta4), path('data/contig/hifi.hap1.txt'), path('data/contig/hifi.hap2.txt'), path('data/contig/ul.hap1.txt'), path('data/contig/ul.hap2.txt'), path('data/contig/hic.hap1.txt'), path('data/contig/hic.hap2.txt')
    tuple val(meta5), path('data/scaffolded/merqury/*'), path('data/scaffolded/stats.hap1.txt'), path('data/scaffolded/stats.hap2.txt'), path('data/scaffolded/busco.hap1.txt'), path('data/scaffolded/busco.hap2.txt'), path('data/scaffolded/snail.hap1.svg'), path('data/scaffolded/snail.hap2.svg'), path('data/scaffolded/blob.hap1.svg'), path('data/scaffolded/blob.hap2.svg'), path('data/scaffolded/pretext.hap1.png'), path('data/scaffolded/pretext.hap2.png')
    tuple val(meta6), path('data/scaffolded/hifi.hap1.txt'), path('data/scaffolded/hifi.hap2.txt'), path('data/scaffolded/ul.hap1.txt'), path('data/scaffolded/ul.hap2.txt'), path('data/scaffolded/hic.hap1.txt'), path('data/scaffolded/hic.hap2.txt')
    tuple val(meta7), path('data/curated/merqury/*'), path('data/contig/curated.hap1.txt'), path('data/contig/curated.hap2.txt'), path('data/curated/busco.hap1.txt'), path('data/curated/busco.hap2.txt'), path('data/curated/snail.hap1.svg'), path('data/curated/snail.hap2.svg'), path('data/curated/blob.hap1.svg'), path('data/curated/blob.hap2.svg'), path('data/curated/pretext.hap1.png'), path('data/curated/pretext.hap2.png')
    tuple val(meta8), path('data/curated/hifi.hap1.txt'), path('data/curated/hifi.hap2.txt'), path('data/curated/ul.hap1.txt'), path('data/curated/ul.hap2.txt'), path('data/curated/hic.hap1.txt'), path('data/curated/hic.hap2.txt')

    output:
    tuple val(meta), path("${prefix}_ear.yaml")                     , emit: yaml
    tuple val(meta), path("${prefix}_ear.pdf")                      , emit: pdf
    tuple val(meta), path("data/", type: 'dir', includeInputs: true), emit: datadir
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_ear, topic: versions

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    def smudgeplot_arg = smudgeplot ? "--smudgeplot $smudgeplot" : ''

    """
    generate_yaml.py \\
        --genomescope $genomescope \\
        $smudgeplot_arg \\
        --contig data/contig/ \\
        --scaffolded data/scaffolded/ \\
        --curated data/curated/ \\
        $args \\
        -o ${prefix}_ear.yaml \\
        $template

    touch ${prefix}_ear.pdf
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ear.yaml
    touch ${prefix}_ear.pdf
    """
}