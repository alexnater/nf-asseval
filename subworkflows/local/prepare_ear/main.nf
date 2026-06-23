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

workflow PREPARE_EAR {

    take:
    ch_evaluation          // channel: [ meta, busco, snail, blob ]
    ch_hic                 // channel: [ meta, snapshot ]
    ch_kmer                // channel: [ meta, genomescope, smudge ]
    ch_merqury             // channel: [ meta, stats, qv, img ]
    ch_depth               // channel: [ meta, depth ]
    ear_yaml               // file: yaml

    main:

    ch_snapshot = ch_hic
        .map { meta, snapshot -> [ [id: meta.ref, sample: meta.sample], snapshot ] }

    ch_depth_stats = ch_depth
        .map { meta, depth -> [ [id: meta.ref, sample: meta.sample], [ meta.type, depth ] ] }
        .groupTuple()
        .map { meta, tuples ->
            def hifi = tuples.find {it[0] == 'hifi'}
            def ul = tuples.find {it[0] == 'ul'}
            def hic = tuples.find {it[0] == 'hic'}
            [ meta, hifi ? hifi[1] : null, ul ? ul[1] : null, hic ? hic[1] : null ]
        }

    ch_kmer_stats = ch_kmer
        .filter { meta, summary, report -> meta.type == 'hifi' }
        .map { meta, summary, report -> [ meta.sample, summary, report ] }

    ch_merqury_stats = ch_merqury
        .map { meta, stats, qv, img -> [ meta.subMap(['sample', 'status']), [ stats, qv, img ].flatten() ] }

    // Add depth output to Busco and branch by haplotype:
    ch_by_type = ch_evaluation
        .map { meta, busco, snail, blob -> [ meta.subMap(['id', 'sample']), meta, busco, snail, blob ] }
        .join(ch_snapshot, failOnDuplicate: true, remainder: true)
        .join(ch_depth_stats, failOnDuplicate: true, remainder: true)
        .branch { key, meta, busco, snail, blob, snapshot, hifi, ul, hic ->
            hap1: meta.type =~ /primary/ || meta.type =~ /hap1/
                return [ meta.subMap(['sample', 'status']), [ busco, snail, blob, snapshot, hifi, ul, hic ] ]  
            hap2: meta.type =~ /alt/ || meta.type =~ /hap2/
                return [ meta.subMap(['sample', 'status']), [ busco, snail, blob, snapshot, hifi, ul, hic ] ]
        }
    
    // Join haplotype pairs with merqury output:
    ch_paired = ch_by_type.hap1
        .join(ch_by_type.hap2, failOnDuplicate: true, remainder: true)
        .join(ch_merqury_stats, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, hap1, hap2, merqury ->
            [ meta + [id: "${meta.sample}_${meta.status}"], [merqury] + [hap1, hap2 ?: [null] * hap1.size()].transpose().flatten() ]
        }

    // Join different assembly stages:
    ch_by_stage = ch_paired
        .map { meta, data -> [ meta.sample, [ meta.status, data.collect { it ?: [] } ] ] }
        .groupTuple()
        .map { sample, tuples ->
            def contig = tuples.find { it[0] == 'contig'}
            def scaffolded = tuples.find { it[0] == 'scaffolded'}
            def curated = tuples.find { it[0] == 'curated'}
            [ sample, contig ? contig[1] : [[]] * 15, scaffolded ? scaffolded[1] : [[]] * 15, curated ? curated[1] : [[]] * 15 ]
        }

    ch_to_ear = ch_by_stage
        .join(ch_kmer_stats, failOnDuplicate: true, remainder: true)
        .multiMap { sample, contig, scaffolded, curated, summary, smudge ->
            yaml:          [ [id: sample], ear_yaml ]
            summary:       [ [id: sample], summary, smudge ]
            contig_stats:  [ [id: "${sample}_contig",     sample: sample, status: 'contig'] ]     + contig[0..8]
            contig_depth:  [ [id: "${sample}_contig",     sample: sample, status: 'contig'] ]     + contig[9..-1]
            scaff_stats:   [ [id: "${sample}_scaffolded", sample: sample, status: 'scaffolded'] ] + scaffolded[0..8]
            scaff_depth:   [ [id: "${sample}_scaffolded", sample: sample, status: 'scaffolded'] ] + scaffolded[9..-1]
            curated_stats: [ [id: "${sample}_curated",    sample: sample, status: 'curated'] ]    + curated[0..8]
            curated_depth: [ [id: "${sample}_curated",    sample: sample, status: 'curated'] ]    + curated[9..-1]
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