//
// Run k-mer analysis
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTK_FASTK                } from '../../../modules/local/fastk/fastk'
include { MERQURYFK_MERQURYFK        } from '../../../modules/nf-core/merquryfk/merquryfk'
include { GENESCOPEFK                } from '../../../modules/nf-core/genescopefk'
include { SMUDGEPLOT                 } from '../../../modules/local/smudgeplot'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_KMER_FK {

    take:
    ch_reads      // channel: [ meta, fastq ]
    ch_fasta_fai  // channel: [ meta, fasta, fai ]
    ch_kmers      // channel: [ meta, db ]
    kmer_size     // value: k-mer size

    main:

    // Group reads by sample and join with pre-existing fastk databases
    ch_to_count = ch_reads
        .map { meta, fastq ->
            def new_meta = [
                id: meta.sample,
                sample: meta.sample,
                type: meta.type,
                samples_per_type: meta.samples_per_type
            ]
            [ meta.sample, new_meta, fastq ]
        }
        .groupTuple(by: [0, 1])
        .join (
            ch_kmers
                .filter { meta, db -> meta.type == 'fastk' && meta.kmer_size == kmer_size }
                .map { meta, db -> [ meta.sample, db ] },
            failOnDuplicate: true,
            remainder: true
        )
        .filter { sample, meta, fastqs, db -> !db }
        .map { sample, meta, fastqs, db -> [ meta, fastqs.flatten() ] }
    
    //
    // MODULE: Run fastk by sample
    //
    FASTK_FASTK (
        ch_to_count,
        kmer_size
    )

    // Join outputs of fastk
    ch_fastk = FASTK_FASTK.out.ktab
        .join(FASTK_FASTK.out.data, failOnDuplicate: true, failOnMismatch: true)
        .join(FASTK_FASTK.out.hist, failOnDuplicate: true, failOnMismatch: true)
        .join(FASTK_FASTK.out.txt, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, ktab, data, hist, txt ->
            [ meta + [kmer_size: kmer_size], ktab, data, hist, txt ]
        }

    //
    // MODULE: Run genescope.fk
    //
    GENESCOPEFK (
        ch_fastk.map { meta, ktab, data, hist, txt -> [ meta, txt ] }
    )

    // Separate assemblies by type:
    ch_by_type = ch_fasta_fai
        .branch { meta, fasta, fai ->
            hap1: meta.type =~ /pri/ || meta.type =~ /hap1/
                return [ meta.subMap(['sample', 'status']), meta, fasta ]  
            hap2: meta.type =~ /alt/ || meta.type =~ /hap2/
                return [ meta.subMap(['sample', 'status']), meta, fasta ]
            other: true
                return [ meta.sample, meta + [hap1: meta.id, hap2: null, status: meta.status], fasta, [] ]
        }

    // Combine haplotypes by sample and status:
    ch_paired = ch_by_type.hap1
        .join(ch_by_type.hap2, failOnDuplicate: true, remainder: true)
        .map { refid, meta, fasta, meta2, fasta2 ->
            [ meta.sample, [hap1: meta.id, hap2: meta2.id, status: meta.status], fasta, fasta2 ]
        }
        .mix(ch_by_type.other)

    // Combine fastk ktabs with reference file
    ch_to_merqury = ch_fastk
        .map { meta, ktab, data, hist, txt -> [ meta.sample, meta, ktab, data, hist, txt ] }
        .combine(ch_paired, by: 0)
        .map { sample, meta, ktab, data, hist, txt, meta2, fasta, fasta2 ->
            [ meta + [hap1: meta2.hap1, hap2: meta2.hap2, status: meta2.status], hist, ktab, data, fasta, fasta2 ]
        }

    //
    // MODULE: Run merqury
    //
    MERQURYFK_MERQURYFK (
        ch_to_merqury,
        [[:], [], []],
        [[:], [], []]
    )

/*
    // Combine fastk ktabs with reference file
    ch_to_merqury = ch_fastk
        .map { meta, ktab, data, hist, txt -> [ meta.sample, meta, ktab, data, hist, txt ] }
        .combine(ch_fasta_fai.map { meta, fasta, fai -> [ meta.sample, meta, fasta, fai ] }, by: 0)
        .map { sample, meta, ktab, data, hist, txt, meta2, fasta, fai ->
            [ meta + [ref: meta2.id], hist, ktab, data, fasta, [] ]
        }

    //
    // MODULE: Run merqury
    //
    MERQURYFK_MERQURYFK (
        ch_to_merqury,
        [[:], [], []],
        [[:], [], []]
    )
*/

    //
    // MODULE: Run smudgeplot
    //
    SMUDGEPLOT (
        ch_fastk.map { meta, ktab, data, hist, txt -> [ meta, ktab, data ] },
        4
    )

    emit:
    summary = GENESCOPEFK.out.summary         // channel: [ meta, summary ]
    stats = MERQURYFK_MERQURYFK.out.stats     // channel: [ meta, stats ]
    qv = MERQURYFK_MERQURYFK.out.qv           // channel: [ meta, qv ]
    images = MERQURYFK_MERQURYFK.out.images   // channel: [ meta, png/pdf ]
    report = SMUDGEPLOT.out.tsv               // channel: [ meta, tsv ]
}