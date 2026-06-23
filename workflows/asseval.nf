/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { WINDOWS_STATS          } from '../modules/local/windows_stats'
include { PLOT_WINDOWS           } from '../modules/local/plot_windows'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { PREPARE_GENOMES        } from '../subworkflows/local/prepare_genomes'
include { RUN_EVALUATION         } from '../subworkflows/local/run_evaluation'
include { RUN_KMER_FK            } from '../subworkflows/local/run_kmer_fk'
include { MAPPABILITY            } from '../subworkflows/local/mappability'
include { MAP_HIC                } from '../subworkflows/local/map_hic'
include { MAP_LONGREADS          } from '../subworkflows/local/map_longreads'
include { MAP_ILLUMINA           } from '../subworkflows/local/map_illumina'
include { BAM_STATS              } from '../subworkflows/local/bam_stats'
include { BAM_DEPTH              } from '../subworkflows/local/bam_depth'
include { RUN_BLOB               } from '../subworkflows/local/run_blob'
include { VARIANT_CALLING        } from '../subworkflows/local/variant_calling'
include { VARIANT_CALLING_GATK   } from '../subworkflows/local/variant_calling_gatk'
include { PREPARE_EAR            } from '../subworkflows/local/prepare_ear'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_asseval_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEVAL {

    take:
    ch_assemblies // channel: samplesheet read in from --assemblies
    ch_reads      // channel: samplesheet read in from --reads
    ch_kmers      // channel: samplesheet read in from --kmers

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    def steps = params.steps ? params.steps.split(',') : []
    if ((steps.contains('stats') || steps.contains('variant_calling')) && !steps.contains('mapping')) {
        steps + ["mapping"]
    }

    //
    // SUBWORKFLOW: prepare_genomes
    // 
    PREPARE_GENOMES ( 
        ch_assemblies.map { meta, fasta, gtf, yaml -> [ meta, fasta, gtf ] },
        ch_reads
    )
    
    ch_fasta_fai = PREPARE_GENOMES.out.genomes
        .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai ] }

    if (steps.contains('qc')) {

        //
        // MODULE: Run FastQC
        //
        FASTQC (
            ch_reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    }

    if (steps.contains('evaluation')) {

        //
        // SUBWORKFLOW: run_evaluation
        //   
        RUN_EVALUATION (
            ch_assemblies.map { meta, fasta, gtf, yaml -> [ meta, fasta ] },
            params.busco_lineage,
            file(params.busco_lineages_path, type: 'dir', checkIfExists: true)
        )
    }

    if (steps.contains('kmer')) {

        //
        // SUBWORKFLOW: run_kmer_fk
        //    
        RUN_KMER_FK (
            ch_reads.filter { meta, fastq -> meta.type == 'hifi' || meta.type == 'illumina' },
            ch_fasta_fai,
            ch_kmers,
            params.kmer_size
        )
    }

    if (steps.contains('mappability')) {

        def kmer_sizes = params.kmer_sizes ? params.kmer_sizes.split(',').collect { it.toInteger() } : []

        //
        // SUBWORKFLOW: mappability
        //
        MAPPABILITY (
            ch_fasta_fai,
            kmer_sizes
        )
        ch_versions = ch_versions.mix(MAPPABILITY.out.versions)
    }

    if (steps.contains('mapping')) {

        // Branch reads by read type:
        ch_reads_bytype = ch_reads
            .branch { meta, reads ->
                longreads: meta.type == 'hifi' || meta.type == 'ont'
                hic: meta.type == 'hic'
                illumina:  true
            }

        //
        // SUBWORKFLOW: map_longreads
        // 
        MAP_LONGREADS (
            ch_reads_bytype.longreads,
            ch_fasta_fai
        )

        //
        // SUBWORKFLOW: map_hic
        //    
        MAP_HIC (
            ch_reads_bytype.hic,
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai, index ] }
        )

        //
        // SUBWORKFLOW: map_illumina
        // 
        MAP_ILLUMINA (
            ch_reads_bytype.illumina,
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai, index ] }
        )

        // Mix all bam files back together
        ch_bam_bai = MAP_LONGREADS.out.bam_bai
            .mix(MAP_HIC.out.bam_bai)
            .mix(MAP_ILLUMINA.out.bam_bai)
    }
    
    if (steps.contains('stats')) {
    
        //
        // SUBWORKFLOW: bam_stats
        //    
        BAM_STATS (
            ch_bam_bai,
            ch_fasta_fai,
            []
        )
        ch_versions = ch_versions.mix(BAM_STATS.out.versions)

        //
        // SUBWORKFLOW: Generate depth per site reports
        //
        BAM_DEPTH (
            ch_bam_bai
        )
        ch_versions = ch_versions.mix(BAM_DEPTH.out.versions)

        //
        // SUBWORKFLOW: run_blob
        //   
        RUN_BLOB (
            ch_assemblies.map { meta, fasta, gtf, yaml -> [ meta, fasta ] },
            RUN_EVALUATION.out.busco_full_table,
            ch_bam_bai.filter { meta, bam, bai -> meta.type == 'hifi' || meta.type == 'illumina' }
        )
    }

    if (steps.contains('variant_calling')) {

        def model_file = params.genotype_model ? file(params.genotype_model, checkIfExists: true) : []
        def config_file = file(params.glnexus_config, checkIfExists: true)

        //
        // SUBWORKFLOW: variant_calling
        //
        VARIANT_CALLING (
            ch_bam_bai.filter { meta, bam, bai -> meta.type == 'hifi' || meta.type == 'ont' },
            ch_fasta_fai,
            [],
            model_file,
            config_file
        )
        ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

        //
        // SUBWORKFLOW: variant_calling_gatk
        //
        VARIANT_CALLING_GATK (
            ch_bam_bai.filter { meta, bam, bai -> meta.type == 'illumina' },
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai, dict ] },
            params.min_contig_length
        )
        ch_versions = ch_versions.mix(VARIANT_CALLING_GATK.out.versions)
    }

    if (steps.contains('report')) {

        // Generate window-wise stats
        ch_summaries = BAM_STATS.out.depth
            .filter { meta, summary -> meta.type == 'hifi' || meta.type == 'ont' }
            .map { meta, summary -> [ groupKey([ref: meta.ref, type: meta.type], meta.samples_per_type), [ meta.sample, summary ] ] }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .map { meta, tuples -> [ meta.target, tuples.collect { it[0] }, tuples.collect { it[1] } ] }

        // Join depth and vcf files per reference
        ch_to_winstats = BAM_DEPTH.out.depth
            .filter { meta, depth, tbi -> meta.type == 'hifi' || meta.type == 'ont' }
            .map { meta, depth, tbi -> [ [ref: meta.ref, type: meta.type], meta, depth, tbi ] }
            .join(VARIANT_CALLING.out.ind_vcf_tbi
                    .filter { meta, vcf, tbi -> (meta.type == 'hifi' || meta.type == 'ont') && meta.caller == 'clair3' }
                    .map { meta, vcf, tbi -> [ [ref: meta.ref, type: meta.type], vcf, tbi ] },
                failOnDuplicate: true,
                failOnMismatch: true
            )
            .join(ch_summaries, failOnDuplicate: true, failOnMismatch: true)
            .map { ref_type, meta, depth, dtbi, vcf, tbi, samples, summaries ->
                [ meta.ref, meta, depth, dtbi, vcf, tbi, samples, summaries ]
            }
            .combine(PREPARE_GENOMES.out.genomes
                    .map { meta, fasta, fai, dict, gtf, index -> [ meta.id, meta, fasta, fai ] },
                by: 0
            )
            .multiMap { ref, meta, depth, dtbi, vcf, tbi, samples, summaries, meta2, fasta, fai ->
                input:     [ [id: 'winstats', type: meta.type, ref: ref, samples: tuple(samples)], depth, dtbi, vcf, tbi, summaries ]
                fasta_fai: [ meta2, fasta, fai ]
            }

        //
        // MODULE: windows_stats
        //
        WINDOWS_STATS (
            ch_to_winstats.input,
            ch_to_winstats.fasta_fai
        )
        ch_versions = ch_versions.mix(WINDOWS_STATS.out.versions.first())

        PLOT_WINDOWS (
            WINDOWS_STATS.out.bed
        )
        ch_versions = ch_versions.mix(PLOT_WINDOWS.out.versions.first())
    }

    if (steps.contains('ear')) {

        def ear_yaml = file(params.ear_yaml, checkIfExists: true)

        ch_evaluation = RUN_EVALUATION.out.busco_summary
            .join(RUN_EVALUATION.out.snail_plot, failOnDuplicate: true, failOnMismatch: true)
            .join(RUN_BLOB.out.blob, failOnDuplicate: true, failOnMismatch: true)

        ch_kmer = RUN_KMER_FK.out.summary
            .join(RUN_KMER_FK.out.report, failOnDuplicate: true, failOnMismatch: true)

        ch_merqury = RUN_KMER_FK.out.stats
            .join(RUN_KMER_FK.out.qv, failOnDuplicate: true, failOnMismatch: true)
            .join(RUN_KMER_FK.out.images, failOnDuplicate: true, failOnMismatch: true)

        //
        // SUBWORKFLOW: generate_ear
        //
        PREPARE_EAR (
            ch_evaluation,
            MAP_HIC.out.snapshot,
            ch_kmer,
            ch_merqury,
            BAM_STATS.out.depth,
            ear_yaml
        )

/*
        ch_depth_stats = BAM_STATS.out.depth
            .map { meta, depth -> [ [id: meta.ref, sample: meta.sample], [ meta.type, depth ] ] }
            .groupTuple()
            .map { meta, tuples ->
                def hifi = tuples.find {it[0] == 'hifi'}
                def ul = tuples.find {it[0] == 'ul'}
                def hic = tuples.find {it[0] == 'hic'}
                [ meta, hifi ? hifi[1] : [], ul ? ul[1] : [], hic ? hic[1] : [] ]
            }

        ch_kmer_stats = RUN_KMER_FK.out.summary
            .join(RUN_KMER_FK.out.report, failOnDuplicate: true, failOnMismatch: true)
            .filter { meta, summary, report -> meta.type == 'hifi' }
            .map { meta, summary, report -> [ meta.sample, summary, report ] }

        ch_merqury = RUN_KMER_FK.out.stats
            .join(RUN_KMER_FK.out.qv, failOnDuplicate: true, failOnMismatch: true)
            .join(RUN_KMER_FK.out.images, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, stats, qv, img -> [ meta.subMap(['sample', 'status']), [ stats, qv, img ].flatten() ] }

        // Add depth output to Busco and branch by haplotype:
        ch_by_type = RUN_EVALUATION.out.busco_summary
            .join(RUN_EVALUATION.out.snail_plot, failOnDuplicate: true, failOnMismatch: true)
            .join(RUN_BLOB.out.blob, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, busco, snail, blob -> [ meta.subMap(['id', 'sample']), meta, busco, snail, blob ] }
            .join(ch_depth_stats, failOnDuplicate: true, remainder: true)
            .branch { key, meta, busco, snail, blob, hifi, ul, hic ->
                hap1: meta.type =~ /primary/ || meta.type =~ /hap1/
                    return [ meta.subMap(['sample', 'status']), busco, snail, blob, hifi, ul, hic ]  
                hap2: meta.type =~ /alt/ || meta.type =~ /hap2/
                    return [ meta.subMap(['sample', 'status']), busco, snail, blob, hifi, ul, hic ]
            }
        
        // Join haplotype pairs with merqury output:
        ch_paired = ch_by_type.hap1
            .join(ch_by_type.hap2, failOnDuplicate: true, remainder: true)
            .join(ch_merqury, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, busco, snail, blob, hifi, ul, hic, busco2, snail2, blob2, hifi2, ul2, hic2, merqury ->
                [ meta + [id: "${meta.sample}_${meta.status}"], [ merqury, busco, busco2, snail, snail2, blob, blob2, hifi, hifi2, ul, ul2, hic, hic2 ] ]
            }

        // Join different assembly stages:
        ch_by_stage = ch_paired
            .map { meta, data -> [ meta.sample, [ meta.status, data ] ] }
            .groupTuple()
            .map { sample, tuples ->
                def contig = tuples.find { it[0] == 'contig'}
                def scaffolded = tuples.find { it[0] == 'scaffolded'}
                def curated = tuples.find { it[0] == 'curated'}
                [ sample, contig ? contig[1] : [], scaffolded ? scaffolded[1] : [], curated ? curated[1] : [] ]
            }

        ch_to_ear = ch_by_stage
            .join(ch_kmer_stats, failOnDuplicate: true, remainder: true)
            .multiMap { sample, contig, scaffolded, curated, summary, smudge ->
                yaml:          [ [id: sample], ear_yaml ]
                summary:       [ [id: sample], summary, smudge ]
                contig_stats:  [ [id: "${sample}_contig",     sample: sample, status: 'contig'] ]     + contig[0..6]
                contig_depth:  [ [id: "${sample}_contig",     sample: sample, status: 'contig'] ]     + contig[7..-1]
                scaff_stats:   [ [id: "${sample}_scaffolded", sample: sample, status: 'scaffolded'] ] + scaffolded[0..6]
                scaff_depth:   [ [id: "${sample}_scaffolded", sample: sample, status: 'scaffolded'] ] + scaffolded[7..-1]
                curated_stats: [ [id: "${sample}_curated",    sample: sample, status: 'curated'] ]    + curated[0..6]
                curated_depth: [ [id: "${sample}_curated",    sample: sample, status: 'curated'] ]    + curated[7..-1]
            }

        //
        // MODULE: generate_ear
        //
        GENERATE_EAR (
            ch_to_ear.yaml,
            ch_to_ear.summary,
            ch_to_ear.contig_stats,
            ch_to_ear.contig_depth,
            ch_to_ear.scaff_stats,
            ch_to_ear.scaff_depth,
            ch_to_ear.curated_stats,
            ch_to_ear.curated_depth
        )
*/

    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'asseval_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
