import sys
import logging
import bisect
import numpy as np
from pathlib import Path
from pysam import VariantFile
from divlib.region import Region
from divlib.depth import Depth

logger = logging.getLogger(__name__)

DTYPE = np.bool


class Hets():
    def __init__(self, region: Region, nsamples: int, init_rows: int=1000):
        self.hets = np.zeros((init_rows, nsamples), dtype=DTYPE)
        self.region = region
        self.pos = np.zeros(init_rows, dtype=np.uint32)
        self.lookup = dict()
        self.npos = 0

    def extend_rows(self, increment: int=1000):
        self.hets = np.append(self.hets, np.zeros((increment, self.hets.shape[1])), axis=0)
        self.pos = np.append(self.pos, np.zeros(increment))

    def shrink_to_size(self):
        self.hets = self.hets[0:self.npos,:]
        self.pos = self.pos[0:self.npos]

    def sort_positions(self):
        self.shrink_to_size()
        sorted_idx = np.argsort(self.pos)
        self.hets = self.hets[sorted_idx,:]
        self.pos = self.pos[sorted_idx]

    def add_row(self, pos: int, hets: np.array=None):
        if self.npos >= self.pos.shape[0]:
            self.extend_rows()
        self.pos[self.npos] = pos
        self.lookup[pos] = self.npos
        if hets is not None:
            if hets.shape[0] != self.hets.shape[1]:
                raise Exception("Length of provided array doesn't match number of samples!")
            self.hets[self.npos,:] = hets
        self.npos += 1

    def set_individual(self, pos: int, sidx: int, het: bool=True):
        if sidx >= self.hets.shape[1]:
            raise Exception("Sample index out of bounds!")
        if not pos in self.lookup:
            self.add_row(pos)
        self.hets[self.lookup[pos], sidx] = het

    def get_window(
            self,
            start_pos: int,
            end_pos: int,
            ) -> np.array:
        wsidx = bisect.bisect_right(self.pos, start_pos)
        weidx = bisect.bisect_right(self.pos, end_pos)
        if wsidx == weidx:
            logger.warning(f"No variant position in window {start_pos}:{end_pos}.")
#        else:
#            logger.info(f"Window {start_pos}:{end_pos} corresponds to position indicies {wsidx} ({int(self.pos[wsidx])}) to {weidx-1} ({int(self.pos[weidx-1])}).")
        return self.hets[wsidx:weidx,:]

    def get_hetsum(self, start_pos: int, end_pos: int, sidx: int) -> np.uint32:
        window = self.get_window(start_pos, end_pos)
        return np.sum(window[:,sidx], axis=0, dtype=np.uint32)

    def get_hetsum_row(self, start_pos: int, end_pos: int) -> np.array:
        window = self.get_window(start_pos, end_pos)
        return np.sum(window, axis=0, dtype=np.uint32)


# Read in vcf file to hets sparse array:
def read_vcf(
        vcf_file: Path,
        depths: Depth,
        region: Region,
        samples: list[str]=None,
        mindepth: int=0,
        mingq: int=0,
        include_indels: bool=True,
        interval: int=1000
        ):
    hets = Hets(region, len(samples))
    with VariantFile(str(vcf_file)) as vcf_in:
        samples = list(vcf_in.header.samples) if samples is None else samples
        if vcf_in.get_tid(region.chrom) < 0:
            return hets, samples
        processed = 0
        for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
            valid_samples = depths.is_valid_row(rec.pos, mindepth)
            gqs = [rec.samples[sample]["GQ"] if rec.samples[sample]["GQ"] is not None else 0 for sample in samples]
            alleles = (al if al is not None else -1 for sample in samples for al in rec.samples[sample]["GT"])
            if rec.rlen == 1 or (rec.rlen > 1 and include_indels):
                hets.add_row(rec.pos, np.array([ True if (is_valid and gq >= mingq and al1 >= 0 and al2 >= 0 and al1 != al2) else False \
                                        for al1, al2, gq, is_valid in zip(alleles, alleles, gqs, valid_samples) ]))
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
        logger.info(f"Extracted {hets.npos} variant sites from VCF file.")
    hets.shrink_to_size()
    return hets, samples

def read_vcf_single(
        vcf_files: list[Path],
        depths: Depth,
        region: Region,
        samples: list[str],
        mindepth: int=0,
        mingq: int=0,
        include_indels: bool=True,
        interval: int=1000
        ):
    sample_dict = {sample: sidx for sidx, sample in enumerate(samples)}
    hets = Hets(region, len(samples))
    for vcf_file in vcf_files:
        with VariantFile(str(vcf_file)) as vcf_in:
            sample = list(vcf_in.header.samples)[0]
            try:
                sidx = sample_dict[sample]
            except:
                raise Exception(f"VCF sample id {sample} is not in list of sample names")
            if vcf_in.get_tid(region.chrom) < 0: continue
            processed = 0
            for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
                is_valid = depths.is_valid(rec.pos, sidx, mindepth)
                gq = rec.samples[sample]["GQ"] if rec.samples[sample]["GQ"] is not None else 0
                alleles = tuple(al if al is not None else -1 for al in rec.samples[sample]["GT"])
                if rec.rlen == 1 or (rec.rlen > 1 and include_indels):
                    is_het = True if (is_valid and gq >= mingq and alleles[0] >= 0 and alleles[1] >= 0 and alleles[0] != alleles[1]) else False
                    hets.set_individual(rec.pos, sidx, is_het)
                processed += 1
                if not processed % interval: logger.info(f"Processed {processed} lines.")
            logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
            logger.info(f"Extracted {hets.npos} variant sites from VCF file.")
    hets.sort_positions()
    return hets, samples
