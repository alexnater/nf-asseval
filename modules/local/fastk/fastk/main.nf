process FASTK_FASTK {
    tag "${meta.id} - ${meta.type}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02c05b2ec421debc83883ef9a211291e3220546c12f7c54cb78e66209cb2797d/data' :
        'community.wave.seqera.io/library/fastk:1.2--4bc70c6cd0d420bd' }"

    input:
    tuple val(meta), path(reads)
    val kvalue

    output:
    tuple val(meta), path("*.hist")                       , emit: hist
    tuple val(meta), path("*.txt")                        , emit: txt
    tuple val(meta), path("*.ktab*")                      , emit: ktab, optional: true
    tuple val(meta), path(".*.ktab*", hidden: true)       , emit: data, optional: true
    tuple val(meta), path("*.{prof,pidx}*", hidden: true) , emit: prof, optional: true
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p tmp

    FastK \\
        -k$kvalue \\
        -T$task.cpus \\
        -M${task.memory.toGiga()} \\
        -N${prefix} \\
        -Ptmp \\
        $args \\
        $reads

    rm -r tmp

    Histex \\
        $args2 \\
        ${prefix}.hist \\
        > ${prefix}.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_ktab = args.contains('-t') ? "touch ${prefix}_fk.ktab .${prefix}_fk.ktab.1" : ''
    def touch_prof = args.contains('-p') ? "touch ${prefix}_fk.prof .${prefix}_fk.pidx.1" : ''
    """
    touch ${prefix}.hist
    touch ${prefix}.txt
    $touch_ktab
    $touch_prof
    """
}
