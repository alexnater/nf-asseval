//
// Map  reads to reference and call variants for Illumina short-read data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                      } from '../../../modules/nf-core/fastp'
include { BWAMEM2_MEM                } from '../../../modules/nf-core/bwamem2/mem'
include { GATK4_MARKDUPLICATES       } from '../../../modules/nf-core/gatk4/markduplicates'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAP_ILLUMINA {

    take:
    ch_reads      // channel: [ meta, fastqs ]
    ch_reference  // channel: [ meta, fasta, fai, index ]

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
    ch_to_map = FASTP.out.reads
        .combine(ch_reference.map { meta, fasta, fai, index -> [ meta, index ] })
        .multiMap { meta, reads, meta2, index ->
            reads: [ meta + [ref: meta2.id], reads ]
            index: [ meta2, index ]
        }

    // Map reads with BWA:
    BWAMEM2_MEM (
        ch_to_map.reads,
        ch_to_map.index,
        [ [:], [] ],
        false
    )

    // group entries by sample:
    ch_bam_bysample = BWAMEM2_MEM.out.bam
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

    // combine bam files with fasta reference and fasta index:
    ch_to_dedup = ch_bam_bysample
        .map { meta, bams -> [ meta.ref, meta.target, bams ] }
        .combine(ch_reference.map { meta, fasta, fai, index ->
            [ meta.id, fasta, fai ]
            },
            by: 0)
        .multiMap { id, meta, bams, fasta, fai ->
            bams:  [ meta, bams ]
            fasta: fasta
            fai:   fai
        }
        
    // merge bam files by sample and mark duplicates:
    GATK4_MARKDUPLICATES (
        ch_to_dedup.bams,
        ch_to_dedup.fasta,
        ch_to_dedup.fai
    )

    // join bam files and the corresponding index files: 
    GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .set { bam_bai }

    emit:
    bam_bai                  // channel: [ val(meta), path(bam), path(bai) ]
}