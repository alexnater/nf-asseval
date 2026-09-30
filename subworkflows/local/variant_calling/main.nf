/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLAIR3                                  } from '../../../modules/local/clair3'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF  } from '../../../modules/nf-core/bcftools/concat'
include { DEEPVARIANT_RUNDEEPVARIANT              } from '../../../modules/nf-core/deepvariant/rundeepvariant'
include { GLNEXUS                                 } from '../../../modules/local/glnexus'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_CALLING {

    take:
    ch_bam_bai     // channel: [ meta, bam, bai ]
    ch_fasta_fai   // channel: [ meta, fasta, fai ]
    bed_file       // BED file with genomic intervals
    model_file     // Clair3 model file
    config_file    // GLnexus config file
    min_length     // int 

    main:

    // Extract list of contigs from each assembly
    ch_regions = ch_fasta_fai
        .map { meta, fasta, fai ->
            def contigs = WorkflowAssEval.getContigs(fai, min_length)
            [ meta + [ncontigs: contigs.size()], contigs.withIndex(1), fasta, fai ]
        }
        .transpose(by: 1)

    // Combine bam files with their reference
    ch_mapped = ch_bam_bai
        .map { meta, bam, bai -> [ meta.ref, meta, bam, bai ] }
        .combine(
            ch_regions.map { meta, contig_with_idx, fasta, fai -> [ meta.id, meta, contig_with_idx, fasta, fai ] },
            by: 0
        )
        .multiMap { ref, meta, bam, bai, meta2, contig_with_idx, fasta, fai ->
            def (contig, idx) = contig_with_idx
            bam_bai:   [ meta + [id: "${meta.id}_${contig}", region: contig, order: idx, ncontigs: meta2.ncontigs], bam, bai, bed_file ]
            fasta_fai: [ meta2, fasta, fai ]
        }

    //
    // MODULE: Run Clair3
    //
    CLAIR3 (
        ch_mapped.bam_bai,
        ch_mapped.fasta_fai,
        ch_mapped.bam_bai.map { meta, bam, bai, bed ->
            [
                meta.type == 'hifi' ? 'hifi' : meta.type == 'ont' ? 'ont' : 'ilmn',
                meta.type == 'hifi' ? 'hifi' : meta.type == 'ont' ? 'ont' : 'ilmn'
            ]
        },
        model_file
    )

    // Concat contig-wise gvcf files
    ch_to_concat = CLAIR3.out.gvcf
        .join(CLAIR3.out.gtbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, gtbi ->
            def new_meta = meta - meta.subMap(['region', 'order', 'ncontigs']) + [id: meta.sample]
            [ groupKey(new_meta, meta.ncontigs), [ meta.order, gvcf, gtbi ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples ->
            def (order, gvcfs, gtbis) = tuples.transpose()
            [ meta, gvcfs, gtbis ]
        }

    //
    // MODULE: Run bcftools concat
    //
    BCFTOOLS_CONCAT_GVCF (
        ch_to_concat,
        false
    )

    ch_from_clair3 = BCFTOOLS_CONCAT_GVCF.out.vcf
        .join(ch_bam_bai, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, bam, bai ->
            def new_meta = [
                id: "joint_${meta.ref}",
                type: meta.type,
                ref: meta.ref,
                caller: 'clair3'
            ]
            [ groupKey(new_meta, meta.samples_per_type), [ meta.sample, gvcf, bam, bai ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples -> 
            def (samples, gvcfs, bams, bais) = tuples.transpose()
            [ meta.target + [samples: tuple(samples)], gvcfs, bams, bais ]
        }

    // Concat contig-wise vcf files
    ch_to_concat2 = CLAIR3.out.vcf
        .join(CLAIR3.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi ->
            [ groupKey(meta - meta.subMap(['region', 'order']), meta.ncontigs) + [id: meta.sample], [ meta.order, vcf, tbi ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples ->
            def (order, vcfs, tbis) = tuples.transpose()
            [ meta, vcfs, tbis ]
        }

    //
    // MODULE: Run bcftools concat
    //
    BCFTOOLS_CONCAT_VCF (
        ch_to_concat2,
        false
    )

    // Join vcf files with index
    BCFTOOLS_CONCAT_VCF.out.vcf
        .join(BCFTOOLS_CONCAT_VCF.out.index, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi ->
            def new_meta = [
                id: "joint_${meta.ref}",
                type: meta.type,
                ref: meta.ref,
                caller: 'clair3'
            ]
            [ groupKey(new_meta, meta.samples_per_type), [ meta.sample, vcf, tbi ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples -> 
            def (samples, vcfs, tbis) = tuples.transpose()
            [ meta.target + [samples: tuple(samples)], vcfs, tbis ]
        }
        .set { ind_vcf_tbi }


/*
    //
    // MODULE: Run DeepVariant
    //
    DEEPVARIANT_RUNDEEPVARIANT (
        ch_mapped.bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, bed_file ]
        },
        ch_mapped.fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        ch_mapped.fasta_fai.map { meta, fasta, fai -> [ meta, fai ] },
        [ [id:'ref'], [] ],
        [ [:], [] ]
    )
        .gvcf
        .join(ch_mapped.bam_bai, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, bam, bai ->
            def new_meta = [
                id: "joint_${meta.ref}",
                type: meta.type,
                ref: meta.ref,
                caller: 'deepvariant'
            ]
            [ groupKey(new_meta, meta.samples_per_type), [ meta.sample, gvcf, bam, bai ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples -> 
            def (samples, gvcfs, bams, bais) = tuples.transpose()
            [ meta.target + [samples: tuple(samples)], gvcfs, bams, bais ]
        }
        .set { ch_from_dv }
*/

    // Mix channels for joint genotyping
    ch_calls = ch_from_clair3
//        .mix(ch_from_dv)
        .multiMap { meta, gvcfs, bams, bais ->
            to_merge: [ meta, gvcfs ]
            bam_bai:  [ meta, bams, bais ]
        }

/*
    //
    // MODULE: Run GLnexus
    //
    GLNEXUS (
        ch_calls.to_merge,
        [ [id: 'regions'], bed_file ],
        ch_calls.to_merge.map { meta, gvcfs -> meta.caller == 'clair3' ? '' : 'DeepVariant' },
        ch_calls.to_merge.map { meta, gvcfs -> meta.caller == 'clair3' ? config_file : [] }
    )
*/

    emit:
//    vcf_tbi  = GLNEXUS.out.vcf_tbi          // channel: [ meta, vcf, tbi ]
    ind_vcf_tbi                             // channel: [ meta, vcf, tbi ]
    bam_bai  = ch_calls.bam_bai             // channel: [ meta, [bam], [bai] ]
}
