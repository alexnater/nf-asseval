//
// Generate report with depth per site
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_DEPTH                       } from '../../../modules/nf-core/samtools/depth'
include { HTSLIB_BGZIPTABIX as BGZIPTABIX_DEPTH } from '../../../modules/nf-core/htslib/bgziptabix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAM_DEPTH {
    take:
    ch_bam_bai        // channel (mandatory): [ val(meta), path(bam), path(bai) ]

    main:

    // Run SAMtools depth over all merged bam files per reference genome:
    ch_to_depth = ch_bam_bai
        .map { meta, bam, bai ->
            def new_meta = [
                id: "depth_${meta.ref}",
                type: meta.type,
                ref: meta.ref
            ]
            [ groupKey(new_meta, meta.samples_per_type), bam ]
        }
        .groupTuple(sort: {a, b -> a.name <=> b.name})
        .map { gkey, bams -> [ gkey.target, bams, [], [] ] }

    SAMTOOLS_DEPTH (
        ch_to_depth
    )

    // bgzip and tabix the depth file:
    BGZIPTABIX_DEPTH (
        SAMTOOLS_DEPTH.out.tsv.map { meta, tsv -> [ meta, tsv, [], [] ] },
        'compress',
        true,
        'tsv'
    )

    BGZIPTABIX_DEPTH.out.output
        .join(BGZIPTABIX_DEPTH.out.index)
        .set { depth }

    emit:
    depth                         // channel: [ val(meta), path(depth.tsv.gz), path(depth.tsv.tbi) ]
}