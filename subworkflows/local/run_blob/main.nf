//
// Run bloobtools
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BLOBTOOLS                  } from '../../../modules/local/blobtools'
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
        .view()
        .multiMap { refid, meta, fasta, table, meta2, bam, bai ->
            fasta: [ meta, fasta, table, [] ]
            bam:   [ meta2, bam ]
        }

    //
    // MODULE: Run blobtools
    //
    BLOBTOOLS (
        ch_to_blob.fasta,
        ch_to_blob.bam
    )

/*
        //
        // MODULE: Run blobtk_create
        //
        BLOBTK_CREATE (
            ch_fasta
                .join(
                    ch_busco,
                    failOnDuplicate:true,
                    failOnMismatch:true
                )
        )

        //
        // MODULE: Run blobtk_depth
        //
        BLOBTK_DEPTH (
            ch_fasta
                .join(
                    ch_busco,
                    failOnDuplicate:true,
                    failOnMismatch:true
                )
        )

        ch_to_plot = ch_fasta
            .join(
                    BLOBTK_CREATE.out.blobdir,
                    failOnDuplicate:true,
                    failOnMismatch:true
                )
            .multiMap { meta, fasta, dir ->
                fasta:   [ meta, fasta ]
                blobdir: dir
            }

        //
        // MODULE: Run blobtk_plot
        //
        BLOBTK_PLOT (
            ch_to_plot.fasta,
            ch_to_plot.blobdir,
            [],
            [],
            'png'
        )

        ch_to_snail = ch_fasta
            .join(
                ch_busco,
                failOnDuplicate:true,
                failOnMismatch:true
            )
            .join(
                    BLOBTK_CREATE.out.blobdir,
                    failOnDuplicate:true,
                    failOnMismatch:true
            )
            .multiMap { meta, fasta, table, dir ->
                fasta:   [ meta, fasta, table ]
                blobdir: [ meta, dir ]
            }

        //
        // MODULE: Run blobtk_snail
        //
        BLOBTK_SNAIL (
            ch_to_snail.fasta,
            ch_to_snail.blobdir,
            'png'
        )
*/
}