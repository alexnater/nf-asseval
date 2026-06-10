process BWAMEM2_MAP_FILTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data'
        : 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113'}"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)
    val(mapq_filter)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('bwamem2'), eval('bwa-mem2 version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_bwamem2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (task.cpus < 4) {
        error "Number of CPUs need to be at least 4"
    }
    
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bwa_cpus = ((task.cpus - 2) / 2) as int

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    mkfifo R1_fifo R2_fifo

    # map forward reads in single-end mode:
    ( bwa-mem2 mem \\
        $args \\
        -t $bwa_cpus \\
        \$INDEX \\
        ${reads[0]} \\
        | filter_five_end.pl > R1_fifo ) &

    # map reverse reads in single-end mode:
    ( bwa-mem2 mem \\
        $args \\
        -t $bwa_cpus \\
        \$INDEX \\
        ${reads[1]} \\
        | filter_five_end.pl > R2_fifo ) &

    # combine filtered single-end bam files:
    ( two_read_bam_combiner_fifo.pl R1_fifo R2_fifo $mapq_filter \\
    | samtools view $args2 --threads 1 -o ${prefix}.bam - ) &

    wait

    rm R1_fifo R2_fifo
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
