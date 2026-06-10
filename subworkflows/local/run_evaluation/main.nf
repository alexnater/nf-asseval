//
// Run assembly evaluation
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GFASTATS                   } from '../../../modules/nf-core/gfastats'
include { QUAST                      } from '../../../modules/nf-core/quast'
include { BUSCO_BUSCO                } from '../../../modules/nf-core/busco/busco'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_EVALUATION {

    take:
    ch_fasta               // channel: [ meta, fasta ]
    busco_lineage          // value
    busco_lineages_dir     // path

    main:

        // Group assemblies by type
        ch_by_type = ch_fasta
            .map { meta, fasta ->
                [ [comp: 'by_type', id: meta.type, type: meta.type], [ meta.id, fasta ] ]
            }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .map { meta, tuples -> [ meta + [labels: tuples.collect { it[0] }], tuples.collect { it[1] } ] }

        // Group assemblies by sample
        ch_by_sample = ch_fasta
            .map { meta, fasta ->
                [ [comp: 'by_sample', id: meta.sample, sample: meta.sample], [ meta.id, fasta ] ]
            }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .map { meta, tuples -> [ meta + [labels: tuples.collect { it[0] }], tuples.collect { it[1] } ] }

        //
        // MODULE: Run gfastats
        //
        GFASTATS (
            ch_fasta,
            [],
            [],
            [],
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ]
        )

        //
        // MODULE: Run quast
        //
        QUAST (
            ch_by_type.mix(ch_by_sample),
            [ [:], [] ],
            [ [:], [] ]
        )

        //
        // MODULE: Run busco
        //
        BUSCO_BUSCO (
            ch_fasta,
            "genome",
            busco_lineage,
            busco_lineages_dir,
            [],
            true
        )

    emit:
    assembly_summary = GFASTATS.out.assembly_summary         // channel: [ meta, summary ]
    busco_summary = BUSCO_BUSCO.out.short_summaries_txt      // channel: [ meta, summary ]
    busco_full_table = BUSCO_BUSCO.out.full_table            // channel: [ meta, table ]
}