#!/usr/bin/env python

import sys
import os
import argparse
import re
from pathlib import Path
import logging
import csv
import yaml

logger = logging.getLogger()


def read_mosdepth(infile: Path) -> dict:
    if not infile.exists():
        return None
    logger.info(f"Processing Mosdepth file {infile} ...")
    with open(infile, 'r', newline='') as inhandle:
        reader = csv.DictReader(inhandle, delimiter='\t')
        for row in reader:
            if row['chrom'] == 'total':
                return float(row['mean'])
    logger.warning(f"No depth information found in file {infile}.")
    return None


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate YAML file for ERGA EAR.",
        epilog="Example: python generate_yaml.py template.yaml -f flagstat.txt -s summary.txt -c coverage.txt -o stats_summary.tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "template",
        type=Path,
        help="YAML template file."
    )
    parser.add_argument(
        "-g",
        "--genomescope",
        type=Path,
        required=True,
        help="GenomeScope summary file"
    )
    parser.add_argument(
        "--smudgeplot",
        type=Path,
        required=False,
        help="Smudgeplot summary file"
    )
    parser.add_argument(
        "--contig",
        type=Path,
        help="Directory with data for contig-level assembly"
    )
    parser.add_argument(
        "--scaffolded",
        type=Path,
        help="Directory with data for scaffolded assembly"
    )
    parser.add_argument(
        "--curated",
        type=Path,
        help="Directory with data for curated assembly"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=Path,
        default="input.yaml",
        help="YAML file for make_EAR.py"
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

    with open(args.template, 'r') as inhandle:
        data = yaml.safe_load(inhandle)
        logger.info(data)

    fn_pattern = re.compile('(hifi|ul|hic)\\.hap(1|2)\\.txt')

    depth_folder = args.contig if args.contig else args.curated
    depths = {}
    for fn in Path(depth_folder).glob('*.txt'):
        if m := re.search(fn_pattern, str(fn)):
            seqtype = m.group(1)
            hidx = int(m.group(2)) - 1
            if not seqtype in depths:
                depths[seqtype] = [0., 0.]
            depths[seqtype][hidx] = read_mosdepth(fn)

    logger.info(f"Depths for assembly: hifi: {depths.get('hifi', 0)}, ONT UL: {depths.get('ul', 0)}, HiC: {depths.get('hic', 0)}.")

    data_list = []
    if 'hifi' in depths:
        data_list.append(f"PacBio HiFi: {depths['hifi']:.2f}x")
    if 'ul' in depths:
        data_list.append(f"ONT UL: {depths['ul']:.2f}x")
    if 'hic' in depths:
        data_list.append(f"HiC: {depths['hic']:.2f}x")
    data['DATA'] = data_list

    if args.contig:
        contig_depth = {}
        for fn in Path(args.contig).glob('*.txt'):
            if m := re.search(fn_pattern, str(fn)):
                seqtype = m.group(1)
                hidx = int(m.group(2)) - 1
                if not seqtype in contig_depth:
                    contig_depth[seqtype] = [0., 0.]
                contig_depth[seqtype][hidx] = read_mosdepth(fn)

        logger.info(f"Depths for contig-level assembly: hifi: {contig_depth.get('hifi', 0)}, ONT UL: {contig_depth.get('ul', 0)}, HiC: {contig_depth.get('hic', 0)}.")

    if args.scaffolded:
        scaffolded_depth = {}
        for fn in Path(args.scaffolded).glob('*.txt'):
            if m := re.search(fn_pattern, str(fn)):
                seqtype = m.group(1)
                hidx = int(m.group(2)) - 1
                if not seqtype in scaffolded_depth:
                    scaffolded_depth[seqtype] = [0., 0.]
                scaffolded_depth[seqtype][hidx] = read_mosdepth(fn)

        logger.info(f"Depths for scaffolded assembly: hifi: {scaffolded_depth.get('hifi', 0)}, ONT UL: {scaffolded_depth.get('ul', 0)}, HiC: {scaffolded_depth.get('hic', 0)}.")

    if args.curated:
        curated_depth = {}
        for fn in Path(args.curated).glob('*.txt'):
            if m := re.search(fn_pattern, str(fn)):
                seqtype = m.group(1)
                hidx = int(m.group(2)) - 1
                if not seqtype in curated_depth:
                    curated_depth[seqtype] = [0., 0.]
                curated_depth[seqtype][hidx] = read_mosdepth(fn)

        logger.info(f"Depths for curated assembly: hifi: {curated_depth.get('hifi', 0)}, ONT UL: {curated_depth.get('ul', 0)}, HiC: {curated_depth.get('hic', 0)}.")

    data['PROFILING']['GenomeScope'] = str(args.genomescope)
    if args.smudgeplot:
        data['PROFILING']['Smudgeplot'] = str(args.smudgeplot)

    with open(args.outfile, 'w') as outhandle:
        print(yaml.safe_dump(data, default_flow_style=False, sort_keys=False), file=outhandle)


if __name__ == "__main__":
    sys.exit(main())
