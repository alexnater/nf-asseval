process CLAIR3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:2.0.3--py311hbc58adc_0' :
        'biocontainers/clair3:2.0.3--py311hbc58adc_0' }"

    input:
    tuple val(meta) , path(bam)  , path(bai), path(bed)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(platform), val(packaged_model)
    path(model)

    output:
    tuple val(meta), path("${prefix}merge_output.vcf.gz"),            emit: vcf
    tuple val(meta), path("${prefix}merge_output.vcf.gz.tbi"),        emit: tbi
    tuple val(meta), path("${prefix}phased_merge_output.vcf.gz"),     emit: phased_vcf, optional: true
    tuple val(meta), path("${prefix}phased_merge_output.vcf.gz.tbi"), emit: phased_tbi, optional: true
    tuple val(meta), path("${prefix}merge_output.gvcf.gz"),           emit: gvcf, optional: true
    tuple val(meta), path("${prefix}merge_output.gvcf.gz.tbi"),       emit: gtbi, optional: true
    tuple val("${task.process}"), val('clair3'), eval('run_clair3.sh --version | sed "s/^Clair3 v//"'), emit: versions_clair3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bed_arg = bed ? "--bed_fn=${bed}" : ''
    def model_path = model ? "${model}" : "/usr/local/bin/models/${packaged_model}"

    """
    run_clair3.sh \\
        --threads=$task.cpus \\
        --sample_name=${meta.id} \\
        --bam_fn=$bam \\
        --ref_fn=$fasta \\
        --output="." \\
        --platform=$platform \\
        --model_path=$model_path \\
        $bed_arg \\
        $args

    # Rename to add prefix
    for file in merge_output.vcf.gz \
            merge_output.vcf.gz.tbi \
            phased_merge_output.vcf.gz \
            phased_merge_output.vcf.gz.tbi \
            merge_output.gvcf.gz \
            merge_output.gvcf.gz.tbi; do
        if [ -e "\$file" ]; then
            mv "\$file" "${prefix}_\${file}"
        fi
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}_phased_merge_output.vcf.gz
    touch ${prefix}_phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}_merge_output.vcf.gz
    touch ${prefix}_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}_merge_output.gvcf.gz
    touch ${prefix}_merge_output.gvcf.gz.tbi
    """
}