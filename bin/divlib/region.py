import sys
import re
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class Region:
    def __init__(self, chrom: str=None, start: int=None, end: int=None, name: str=None):
        if chrom is None or start is None or end is None:
            self.defined = False
        else:
            self.defined = True
        if self.defined and start > end:
            raise ValueError (f"Start value {start} is larger than end value {end}!")
        self.chrom = chrom
        self.startpos = start
        self.endpos = end
        self.name = name if name else f"{chrom}:{start}-{end}"

    @classmethod
    def from_string(cls, region_str: str):
        tmp = region_str.split(':')
        if len(tmp) != 2:
            raise ValueError (f"Invalid region string {region_str}!")
        positions = tuple(int(x) for x in tmp[1].split('-'))
        if len(positions) != 2:
            raise ValueError (f"Invalid region string {region_str}!")
        return cls(tmp[0], positions[0], positions[1])

    def is_in_interval(self, chrom: str, pos: int) -> bool:
        if not chrom == self.chrom:
            return False
        if pos >= self.startpos and pos <= self.endpos:
            return True
        return False
    
    def __str__(self):
        return f"{self.chrom}_{self.start}-{self.end}"
    
    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    @property    
    def start(self) -> int:
        if self.startpos is None:
            return 1
        else:
            return self.startpos
        
    @property    
    def end(self) -> int:
        return self.endpos

    @property    
    def length(self) -> int:
        if self.defined:
            return self.endpos - self.startpos + 1
        else:
            return None


def get_regions(bed_file: Path=None, gtf_file: Path=None, fai_file: Path=None, region_strings: "list[str]"=None) -> "list[Region]":
    regions = []
    if region_strings:
        for region_str in region_strings:
            try:
                region = Region.from_string(region_str)
            except ValueError as e:
                logger.critical(e)
                sys.exit(1)
            regions.append(region)
    if bed_file:
        regions.extend(read_bed(bed_file))
    if gtf_file:
        regions.extend(read_gtf(gtf_file))
    if fai_file:
        regions.extend(read_fai(fai_file))
    if len(regions) == 0:
        logger.error(f"No valid region provided!")
        sys.exit(1)
    return regions

def read_fai(fai_file: Path) -> "list[Region]":
    regions = []
    with open(fai_file, 'r') as fai_in:
        for line in fai_in:
            fields = line.strip().split('\t')
            if not len(fields) == 5:
                raise Exception("Invalid format of line in fasta index!")
            regions.append(Region(fields[0], 1, int(fields[1])))
    return regions

def read_bed(bed_file: Path, interval: int=1000) -> "list[Region]":
    regions = []
    with open(bed_file, 'r') as inhandle:
        processed = 0
        for line in inhandle:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                logger.error(f"Wrongly formatted line in BED file!")
                sys.exit(1)
            regions.append(Region(fields[0], int(fields[1]) + 1, int(fields[2])))
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of BED file.")
    return regions

def read_gtf(gtf_file: Path, interval: int=1000) -> "list[Region]":
    regions = []
    with open(gtf_file, 'r') as inhandle:
        idpattern = re.compile(r'gene_id "(.+)"')
        processed = 0
        for line in inhandle:
            fields = line.strip().split('\t')
            if len(fields) != 9:
                logger.error(f"Wrongly formatted line in GTF file!")
                sys.exit(1)
            m = idpattern.match(fields[8])
            try:
                region_name = m.group(1)
            except AttributeError:
                logger.error(f"Feature doesn't have a valid tag ({fields[8]}). Skipping line.")
                continue
            regions.append(Region(fields[0], int(fields[3]), int(fields[4]), region_name))
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of GTF file.")
    return regions
