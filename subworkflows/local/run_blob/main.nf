//
// Run bloobtools
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BLOBTOOLS_BLOB             } from '../../../modules/local/blobtools/blob'
/*
include { BLOBTK_CREATE              } from '../../../modules/nf-core/blobtk/create'
include { BLOBTK_DEPTH               } from '../../../modules/nf-core/blobtk/depth'
include { BLOBTK_PLOT                } from '../../../modules/nf-core/blobtk/plot'
include { BLOBTK_SNAIL               } from '../../../modules/nf-core/blobtk/snail'
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_BLOB {

    take:
    ch_fasta               // channel: [ meta, fasta ]
    ch_busco               // channel: [ meta, table ]
    ch_bam_bai             // channel: [ meta, bam, bai ]

    main:

    // Combine references with bam files and busco tables
    ch_to_blob = ch_fasta
        .join(ch_busco, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, fasta, table -> [ [id: meta.id, sample: meta.sample], meta, fasta, table ] }
        .combine(ch_bam_bai.map { meta, bam, bai -> [ [id: meta.ref, sample: meta.sample], meta, bam, bai ] }, by: 0)
        .multiMap { refid, meta, fasta, table, meta2, bam, bai ->
            fasta: [ meta, fasta, table, [] ]
            bam:   [ meta2, bam ]
        }

    //
    // MODULE: Run blobtools blob
    //
    BLOBTOOLS_BLOB (
        ch_to_blob.fasta,
        ch_to_blob.bam
    )

    emit:
    blob = BLOBTOOLS_BLOB.out.blob        // channel: [ meta, png/pdf ]
}