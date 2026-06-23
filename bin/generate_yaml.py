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

    data['DATA'] = []
    if 'hifi' in depths:
        data['DATA'].append(f"PacBio HiFi: {depths['hifi']:.2f}x")
    if 'ul' in depths:
        data['DATA'].append(f"ONT UL: {depths['ul']:.2f}x")
    if 'hic' in depths:
        data['DATA'].append(f"HiC: {depths['hic']:.2f}x")

    data['PROFILING']['GenomeScope'] = str(args.genomescope)
    if args.smudgeplot:
        data['PROFILING']['Smudgeplot'] = str(args.smudgeplot)

    if args.scaffolded:
        hap1 = data['ASSEMBLIES']['Pre-curation']['hap1']
        hap1['gfastats--nstar-report_txt'] = Path(args.scaffolded, 'stats.hap1.txt')
        hap1['busco_short_summary_txt'] = Path(args.scaffolded, 'busco.hap1.txt')
        hap1['merqury_folder'] = Path(args.scaffolded, 'merqury')

        hap2 = data['ASSEMBLIES']['Pre-curation']['hap2']
        hap2['gfastats--nstar-report_txt'] = Path(args.scaffolded, 'stats.hap2.txt')
        hap2['busco_short_summary_txt'] = Path(args.scaffolded, 'busco.hap2.txt')
        hap2['merqury_folder'] = Path(args.scaffolded, 'merqury')

    if args.curated:
        hap1 = data['ASSEMBLIES']['Curated']['hap1']
        hap1['gfastats--nstar-report_txt'] = Path(args.curated, 'stats.hap1.txt')
        hap1['busco_short_summary_txt'] = Path(args.curated, 'busco.hap1.txt')
        hap1['merqury_folder'] = Path(args.curated, 'merqury')
        hap1['hic_FullMap_png'] = Path(args.curated, 'hic.hap1.png')
        hap1['blobplot_cont_png'] = Path(args.curated, 'blob.hap1.svg')

        hap2 = data['ASSEMBLIES']['Curated']['hap2']
        hap2['gfastats--nstar-report_txt'] = Path(args.curated, 'stats.hap2.txt')
        hap2['busco_short_summary_txt'] = Path(args.curated, 'busco.hap2.txt')
        hap2['merqury_folder'] = Path(args.curated, 'merqury')
        hap2['hic_FullMap_png'] = Path(args.curated, 'hic.hap2.png')
        hap2['blobplot_cont_png'] = Path(args.curated, 'blob.hap2.svg')


    with open(args.outfile, 'w') as outhandle:
        print(yaml.safe_dump(data, default_flow_style=False, sort_keys=False), file=outhandle)


if __name__ == "__main__":
    sys.exit(main())
