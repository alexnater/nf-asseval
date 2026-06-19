//
// Generate ERGA Assembly Report
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENERATE_EAR               } from '../../../modules/local/generate_ear'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_EAR {

    take:
    ch_busco               // channel: [ meta, summary ]
    ch_blob                // channel: [ meta, blob ]
    ch_summary             // channel: [ meta, summary ]
    ch_report              // channel: [ meta, report ]
    ch_stats               // channel: [ meta, stats ]
    ch_qv                  // channel: [ meta, qv ]
    ch_images              // channel: [ meta, images ]
    ch_depth               // channel: [ meta, depth ]
    ear_yaml               // file: yaml

    main:

    ch_depth_stats = ch_depth
        .map { meta, depth -> [ [id: meta.ref, sample: meta.sample], [ meta.type, depth ] ] }
        .groupTuple()
        .map { meta, tuples ->
            def hifi = tuples.find {it[0] == 'hifi'}
            def ul = tuples.find {it[0] == 'ul'}
            def hic = tuples.find {it[0] == 'hic'}
            [ meta, hifi ? hifi[1] : [], ul ? ul[1] : [], hic ? hic[1] : [] ]
        }

    ch_kmer_stats = ch_summary
        .join(ch_report, failOnDuplicate: true, failOnMismatch: true)
        .filter { meta, summary, report -> meta.type == 'hifi' }
        .map { meta, summary, report -> [ meta.sample, summary, report ] }

    ch_merqury = ch_stats
        .join(ch_qv, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_images, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, stats, qv, img -> [ meta.subMap(['sample', 'status']), [ stats, qv, img ].flatten() ] }

    // Add depth output to Busco and branch by haplotype:
    ch_by_type = RUN_EVALUATION.out.busco_summary
        .join(ch_blob, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, busco, blob -> [ meta.subMap(['id', 'sample']), meta, busco, blob ] }
        .join(ch_depth_stats, failOnDuplicate: true, remainder: true)
        .branch { key, meta, busco, blob, hifi, ul, hic ->
            hap1: meta.type =~ /primary/ || meta.type =~ /hap1/
                return [ meta.subMap(['sample', 'status']), busco, blob, hifi, ul, hic ]  
            hap2: meta.type =~ /alt/ || meta.type =~ /hap2/
                return [ meta.subMap(['sample', 'status']), busco, blob, hifi, ul, hic ]
        }
    
    // Join haplotype pairs with merqury output:
    ch_paired = ch_by_type.hap1
        .join(ch_by_type.hap2, failOnDuplicate: true, remainder: true)
        .join(ch_merqury, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, busco, blob, hifi, ul, hic, busco2, blob2, hifi2, ul2, hic2, merqury ->
            [ meta + [id: "${meta.sample}_${meta.status}"], [ merqury, busco, busco2, blob, blob2, hifi, hifi2, ul, ul2, hic, hic2 ] ]
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
            contig_stats:  [ [id: "${sample}_contig",     sample: sample, status: 'contig'] ]     + contig[0..4]
            contig_depth:  [ [id: "${sample}_contig",     sample: sample, status: 'contig'] ]     + contig[5..-1]
            scaff_stats:   [ [id: "${sample}_scaffolded", sample: sample, status: 'scaffolded'] ] + scaffolded[0..4]
            scaff_depth:   [ [id: "${sample}_scaffolded", sample: sample, status: 'scaffolded'] ] + scaffolded[5..-1]
            curated_stats: [ [id: "${sample}_curated",    sample: sample, status: 'curated'] ]    + curated[0..4]
            curated_depth: [ [id: "${sample}_curated",    sample: sample, status: 'curated'] ]    + curated[5..-1]
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
}