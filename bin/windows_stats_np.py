#!/usr/bin/env python

import sys
import argparse
import logging
import csv
import numpy as np
from pathlib import Path
from collections import Counter
from pysam import VariantFile
from divlib.region import Region, get_regions
from divlib.depth import Depth, read_depth, get_samples_from_depth_file
from divlib.heterozygosity import Hets, read_vcf, read_vcf_single
from divlib.fasta import get_sequence

logger = logging.getLogger()


def calculate_window_stats(
        region: Region,
        depths: Depth,
        hets: Hets,
        refseq: str,
        windowsize: int,
        stepsize: int,
        mindepth: int,
        total_dp: list[float]
        ):
    cgcs, cats = 0., 0.
    winvalues = []
    winstart = region.start - 1
    while winstart < region.end:
        winend = winstart + windowsize
        winend = winend if winend < region.end else region.end
        startpos = winstart + 1
        gc, gc_skew, at_skew = float('nan'), float('nan'), float('nan')
        if not refseq is None:
            seqstart = startpos - region.start
            seqend = winend - region.start + 1
            counter = Counter(refseq[seqstart:seqend].upper())
            totvalid = (counter['A'] + counter['T'] + counter['G'] + counter['C'])
            if totvalid > 0:
                gc = (counter['G'] + counter['C']) / totvalid
            if (counter['G'] + counter['C']) > 0:
                gc_skew = (counter['G'] - counter['C']) / (counter['G'] + counter['C'])
                cgcs += gc_skew
            if (counter['A'] + counter['T']) > 0:
                at_skew = (counter['A'] - counter['T']) / (counter['A'] + counter['T'])
                cats += at_skew
        logger.info(f"Working on window {region.chrom}:{startpos}-{winend} ...")
        nvalid = depths.get_nvalid_row(startpos, winend, mindepth)
        covsums = depths.get_colsum_row(startpos, winend)
        hetsums = hets.get_hetsum_row(startpos, winend)
        winvalues.append([region.chrom, winstart, winend, gc, gc_skew, cgcs, at_skew, cats, \
                          nvalid, hetsums, \
                          covsums / (winend - winstart), \
                          covsums / (np.array(total_dp, dtype=np.float64) * (winend - winstart))])
        winstart += stepsize
    return winvalues


def read_vcf_full(
        vcf_file: Path,
        depths: Depth,
        region: Region,
        samples: list[str]=None,
        mindepth: int=0,
        include_indels: bool=True,
        interval: int=1000
        ):
    hets = np.zeros((region.length, len(samples)), dtype=np.bool)
    with VariantFile(str(vcf_file)) as vcf_in:
        samples = list(vcf_in.header.samples) if samples is None else samples
        if vcf_in.get_tid(region.chrom) < 0:
            return samples, hets
        processed = 0
        for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
            valid_samples = [depth >= mindepth for depth in depths.get_position_row(rec.pos)]
            alleles = (al if al is not None else -1 for sample in samples for al in rec.samples[sample]["GT"])
            if rec.rlen == 1 or (rec.rlen > 1 and include_indels):
                hets[rec.pos] = np.array([ True if (is_valid and al1 >= 0 and al2 >= 0 and al1 != al2) else False \
                                        for al1, al2, is_valid in zip(alleles, alleles, valid_samples) ])
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
    return samples, hets


def read_mosdepth(infiles: list[Path]) -> tuple[dict, list, dict]:
    total = [0.] * len(infiles)
    depth = {}
    chrlen = {}
    for fidx, infile in enumerate(infiles):
        with open(infile, 'r', newline='') as inhandle:
            reader = csv.DictReader(inhandle, delimiter='\t')
            for row in reader:
                if row['chrom'] == 'total':
                    total[fidx] = float(row['mean'])
                if not row['chrom'] in chrlen:
                    chrlen[row['chrom']] = int(row['length'])
                elif row['chrom'] != 'total' and int(row['length']) != chrlen[row['chrom']]:
                    logger.error(f"ERROR: Unequal length of reference contig {row['chrom']}: {row['length']} vs. {chrlen[row['chrom']]}")
                if not row['chrom'] in depth:
                    depth[row['chrom']] = [0.] * len(infiles)
                depth[row['chrom']][fidx] = float(row['mean'])
    return depth, total, chrlen


def print_summary(outfile: Path, depth: dict, total: list[float], chrom_lengths: dict, samples: list[str]):
    with open(outfile, 'w') as outhandle:
        print("contig", "length", "\t".join(f"raw_{sample.split('_')[0]}" for sample in samples), \
              "\t".join(f"corr_{sample.split('_')[0]}" for sample in samples), "F/M", sep="\t", file=outhandle)
        for contig, dps in depth.items():
            print(contig, chrom_lengths[contig], "\t".join(f"{dp:.2f}" for dp in dps), \
                  "\t".join(f"{dp / tot:.2f}" for dp, tot in zip(dps, total)), \
                    f"{(dps[2] / dps[1]) * (total[1] / total[2]):.2f}" if len(samples) >= 3 and dps[1] > 0 else float("nan"), \
                    sep="\t", file=outhandle)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculates windowwise corrected mean depth per sample by region.",
        epilog="Example: python windows_stats.py --depth depth.tsv --vcf variants.vcf -o out1 --samples sample1 sample2 sample3 --bed regions.bed --windowsize 10000 --stepsize 2000",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--depth",
        type=Path,
        required=True,
        help="Site-wise coverage information in bgzipped and tabixed SAMtools depth format."
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        nargs='+',
        required=True,
        help="VCF file with joint variant calls or list of single-sample VCF files, bgzipped and tabixed."
    )
    parser.add_argument(
        "-o",
        "--outprefix",
        type=str,
        required=True,
        help="Prefix for output files."
    )
    parser.add_argument(
        "-s",
        "--samples",
        nargs='+',
        type=str,
        required=False,
        help="List of sample IDs."
    )
    parser.add_argument(
        "-r",
        "--regions",
        nargs='+',
        type=str,
        required=False,
        help="List of regions in chrom:start-end format to work on.",
    )
    parser.add_argument(
        "--summaries",
        nargs='+',
        type=Path,
        required=False,
        help="List of mosdepth summary files.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=False,
        help="Reference genome FASTA file."
    )
    parser.add_argument(
        "--fai",
        type=Path,
        required=False,
        help="Fasta index of reference genome."
    )
    parser.add_argument(
        "--bed",
        type=Path,
        required=False,
        help="BED file with regions to process."
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        required=False,
        help="GTF file with regions to process."
    )
    parser.add_argument(
        "--minlength",
        type=int,
        help="Minium length of region to consider.",
        default=1000000
    )
    parser.add_argument(
        "-d",
        "--mindepth",
        type=int,
        help="Minium depth needed to be considered as valid genotype call.",
        default=5
    )
    parser.add_argument(
        "-q",
        "--mingq",
        type=int,
        help="Minium genotype quality needed to be considered as valid genotype call.",
        default=0
    )
    parser.add_argument(
        "--include_indels",
        help="Include INDEL variants.",
        action='store_true'
    )
    parser.add_argument(
        "-w",
        "--windowsize",
        type=int,
        help="Size of windows.",
        default=100
    )
    parser.add_argument(
        "--stepsize",
        type=int,
        help="Stepsize for sliding windows.",
        default=100
    )
    parser.add_argument(
        "--total_depth",
        nargs='+',
        type=float,
        required=False,
        help="Total genomic mean depth for each sample.",
    )
    parser.add_argument(
        "-i",
        "--interval",
        type=int,
        help="The desired interval of progress updates.",
        default=1000
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING"
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    regions = [region for region in get_regions(args.bed, args.gtf, args.fai, args.regions) if region.length >= args.minlength]
    samples = args.samples if args.samples else get_samples_from_depth_file(args.depth)
    total_depth = args.total_depth if args.total_depth else [1.] * len(samples)

    if args.summaries:
        chrom_depth, mean_dp, chrom_lengths = read_mosdepth(args.summaries)
        print_summary(Path(f"{args.outprefix}.summary.tsv"), chrom_depth, mean_dp, chrom_lengths, samples)
        if not args.total_depth:
            logger.info(f"Setting mean depth to {mean_dp}.")
            total_depth = mean_dp

    with open(Path(f"{args.outprefix}_win{args.windowsize}_step{args.stepsize}.bed"), 'w') as outhandle:
        print("#chrom", "start", "end", "gc_content", "gc_skew", "cgc_skew", "at_skew", "cat_skew", \
            "\t".join(f"{sample.split('_')[0]}_val" for sample in samples), \
            "\t".join(f"{sample.split('_')[0]}_het" for sample in samples), \
            "\t".join(f"{sample.split('_')[0]}_raw" for sample in samples), \
            "\t".join(f"{sample.split('_')[0]}_corr" for sample in samples), sep="\t", file=outhandle)

        for region in regions:
            if args.fasta:
                try:
                    refseq = get_sequence(args.fasta, region)
                except KeyError as e:
                    logger.critical(f"Couldn't find {region} in reference sequence!")
                    sys.exit(1)
            else:
                refseq = None
            
            logger.info(f"Working on region {region} ...")
            depths = read_depth(args.depth, region, samples, args.interval)
            if len(args.vcf) == 1:
                hets, _ = read_vcf(args.vcf[0], depths, region, samples, args.mindepth, args.mingq, args.include_indels, args.interval)
            else:
                hets, _ = read_vcf_single(args.vcf, depths, region, samples, args.mindepth, args.mingq, args.include_indels, args.interval)
            winvalues = calculate_window_stats(region, depths, hets, refseq, args.windowsize, \
                                                args.stepsize, args.mindepth, total_depth)
            for chrom, start, end, gc, gc_skew, cgcs, at_skew, cats, valid, hets, dps_raw, dps_corr in winvalues:
                print(chrom, str(start), str(end), f"{gc:.2f}", f"{gc_skew:.2f}", f"{cgcs:.2f}", f"{at_skew:.2f}", f"{cats:.2f}", \
                    "\t".join(f"{val}" for val in valid), "\t".join(f"{het}" for het in hets), \
                    "\t".join(f"{cov:.4f}" for cov in dps_raw), \
                    "\t".join(f"{cov:.4f}" for cov in dps_corr), sep="\t", file=outhandle)


if __name__ == "__main__":
    sys.exit(main())
