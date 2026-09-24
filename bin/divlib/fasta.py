import sys
import logging
from pathlib import Path
from pysam import FastaFile
from divlib.region import Region

logger = logging.getLogger(__name__)


def get_sequence(fasta_file: Path, region: Region) -> str:
    with FastaFile(str(fasta_file)) as fasta_in:
        seq = fasta_in.fetch(region.chrom, region.start-1, region.end)
    return seq
