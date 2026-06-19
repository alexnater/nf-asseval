#!/usr/bin/env python

import sys
import os
import argparse
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

    if args.contig:
        hifi1 = read_mosdepth(Path(args.contig, 'hifi.hap1.txt'))
        hifi2 = read_mosdepth(Path(args.contig, 'hifi.hap2.txt'))
        ul1 = read_mosdepth(Path(args.contig, 'ul.hap1.txt'))
        ul2 = read_mosdepth(Path(args.contig, 'ul.hap2.txt'))
        hic1 = read_mosdepth(Path(args.contig, 'hic.hap1.txt'))
        hic2 = read_mosdepth(Path(args.contig, 'hic.hap2.txt'))
        logger.info(f"Depth for contig-level assembly: {hifi1}, {hifi2}, {ul1}, {ul2}, {hic1}, {hic2}.")

    with open(args.outfile, 'w') as outhandle:
       yaml.safe_dump(data)


if __name__ == "__main__":
    sys.exit(main())
