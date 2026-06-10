//
// Map Hi-C reads to reference and produce PretextView file
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                      } from '../../../modules/nf-core/fastp'
include { BWAMEM2_MAP_FILTER         } from '../../../modules/local/bwamem2/map_filter'
include { GATK4_MARKDUPLICATES       } from '../../../modules/nf-core/gatk4/markduplicates'
include { PRETEXTMAP                 } from '../../../modules/nf-core/pretextmap'
include { PRETEXTSNAPSHOT            } from '../../../modules/nf-core/pretextsnapshot'
include { YAHS                       } from '../../../modules/nf-core/yahs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAP_HIC {

    take:
    ch_reads        // channel: [ meta, fastqs ]
    ch_reference    // channel: [ meta, fasta, fai, index ]

    main:

    //
    // MODULE: Run fastp
    //
    FASTP (
        ch_reads.map { meta, fastqs -> [ meta, fastqs, [] ] },
        false,
        false,
        false
    )

    // Prepare channel that joins reads with reference inidices:
    FASTP.out.reads
        .combine(ch_reference.map { meta, fasta, fai, index -> [ meta, index ] })
        .multiMap { meta, reads, meta2, index ->
            reads: [ meta + [ref: meta2.id], reads ]
            index: [ meta2, index ]
        }
        .set { ch_to_map }

    // Map reads with BWA:
    BWAMEM2_MAP_FILTER (
        ch_to_map.reads,
        ch_to_map.index,
        params.mapq_filter
    )

    // group entries by sample:
    BWAMEM2_MAP_FILTER.out.bam
        .map { meta, bam ->
            def new_meta = [
                id: meta.sample,
                sample: meta.sample,
                type: meta.type,
                ref: meta.ref,
                samples_per_type: meta.samples_per_type
            ]
            [ groupKey(new_meta, meta.runs_per_sample), bam ]
        }
        .groupTuple(sort: true)
        .set { ch_bam_bysample }

    // combine bam files with fasta reference and fasta index:
    ch_bam_bysample
        .map { meta, bams -> [ meta.ref, meta.target, bams ] }
        .combine(ch_reference.map { meta, fasta, fai, index ->
            [ meta.id, fasta, fai ]
            },
            by: 0)
        .multiMap { id, meta, bams, fasta, fai ->
            bams:  [ meta, bams ]
            fasta: fasta
            fai:   fai
        }.set { ch_to_dedup }
        
    // merge bam files by sample and mark duplicates:
    GATK4_MARKDUPLICATES (
        ch_to_dedup.bams,
        ch_to_dedup.fasta,
        ch_to_dedup.fai
    )

    // join bam files and the corresponding index files: 
    GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_bam_bai }

/*  THIS DOESN'T WORK
    // join all hic reads:
    ch_bam_bai
        .map { meta, bam, bai -> [ [id: meta.ref], bam, bai ] }
        .groupTuple()
        .set { ch_bams_by_ref }
*/

    //
    // MODULE: Run pretextmap
    //
    PRETEXTMAP (
        ch_bam_bai.map { meta, bam, bai -> [ meta, bam ] },
        [ [:], [], [] ]
    )

    //
    // MODULE: Run pretextsnapshot
    //
    PRETEXTSNAPSHOT (
        PRETEXTMAP.out.pretext
            .map { meta, pretext -> [ meta, pretext, [] ] }
    )

    // combine bam files with fasta reference and fasta index:
    ch_to_yahs = ch_bam_bai
        .map { meta, bam, bai -> [ meta.ref, meta, bam, bai ] }
        .combine(
            ch_reference.map { meta, fasta, fai, index ->
                [ meta.id, fasta, fai ]
            },
            by: 0
        )
        .map { id, meta, bam, bai, fasta, fai ->
            [ meta, fasta, fai, bam, [] ]
        }

    //
    // MODULE: Run yahs
    //
    YAHS (
        ch_to_yahs
    )

    emit:
    bam_bai  = ch_bam_bai                 // channel: [ val(meta), path(bam), path(bai) ]
    pretext  = PRETEXTMAP.out.pretext     // channel: [ val(meta), path(pretext) ]
}