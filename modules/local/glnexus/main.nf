process GLNEXUS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c385aa04ad92c131a0cf00d4070a94aafb24eca9ce323ce49f91b74ee557b58f/data' :
        'community.wave.seqera.io/library/bcftools_glnexus_htslib:626d3efb138422e9' }"

    input:
    tuple val(meta), path(gvcfs)
    tuple val(meta2), path(bed)
    val(preset)
    path(config)

    output:
    tuple val(meta), path("*.bcf")                  , emit: bcf
    tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf_tbi
    tuple val("${task.process}"), val('glnexus'), eval("glnexus_cli 2>&1 | head -n 1 | sed 's/^.*release v//; s/ .*\$//'"), emit: versions_glnexus, topic: versions
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'"), emit: versions_bcftools, topic: versions
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | sed 's/^.*Version: //; s/ .*\$//'"), emit: versions_tabix, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args2 ?: ''
    def args4 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--bed ${bed}" : ""
    def config_command = preset ? "--config $preset" : config ? "--config $config" : ""

    // Make list of GVCFs to merge
    def input = gvcfs.collect { it.toString() }
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Glnexus] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    glnexus_cli \\
        --threads $task.cpus \\
        --mem-gbytes $avail_mem \\
        $regions \\
        $config_command \\
        $args \\
        ${input.join(' ')} \\
        > ${prefix}.bcf

    bcftools \\
        view \\
        $args2 \\
        ${prefix}.bcf \\
        | bgzip \\
            --threads ${task.cpus} \\
            $args3 \\
            > ${prefix}.vcf.gz

    tabix \\
        --threads ${task.cpus} \\
        -p vcf \\
        $args4 \\
        ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf
    """
}
