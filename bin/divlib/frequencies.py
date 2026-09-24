import sys
import logging
from pathlib import Path
import numpy as np
from pysam import TabixFile, asTuple
from divlib.general import Base
from divlib.region import Region

logger = logging.getLogger(__name__)

DTYPE = np.float64


class Frequencies():
    def __init__(self, region: Region, pops: list[str]):
        self.freqs = np.empty((region.length, len(pops), 4), dtype=DTYPE)
        self.freqs.fill(np.nan)
        self.nind = np.zeros((region.length, len(pops)), dtype=np.uint8)
        self.region = region
        self.pops = pops

    def is_valid(self, pos: int, pidx: int) -> bool:
        return np.all(np.isnan(self.freqs[pos-self.region.start, pidx,:]))

    def get_position(self, pos: int, pidx: int) -> float:
        return self.freqs[pos-self.region.start, pidx,:]
    
    def get_position_row(self, pos: int) -> np.array:
        return self.freqs[pos-self.region.start,:,:]

    def set_position(self, pos: int, pidx: int, freqs: list[float]):
        self.freqs[pos-self.region.start, pidx,:] = freqs

    def get_pop(self, spos: int, epos: int, pidx: int) -> np.array:
        sidx = spos-self.region.start
        eidx = epos-self.region.start + 1
        return self.freqs[sidx:eidx,pidx,:]
    
    def get_diversity(
            self,
            spos: int,
            epos: int,
            minprop: float=0.1
            ) -> tuple[float,float,float,float,float,int]:
        sidx = spos-self.region.start
        eidx = epos-self.region.start + 1
        window = self.freqs[sidx:eidx,:,:]
        nvalid = np.sum(~np.isnan(window[:,:,0]), 0)
        nvalid_tot = np.sum(~(np.isnan(window[:,0,0]) | np.isnan(window[:,1,0])))
        minvalid = (eidx - sidx) * minprop
        if nvalid_tot < minvalid:
            return np.nan, np.nan, np.nan, np.nan, np.nan, nvalid_tot
        spi1 = 1 - np.sum(np.square(window[:,0,:]), axis=1)
        spi2 = 1 - np.sum(np.square(window[:,1,:]), axis=1)
        spitot = 1 - np.sum(np.square((window[:,0,:] + window[:,1,:])/2), axis=1)
        sdxy = 1 - np.sum(window[:,0,:] * window[:,1,:], axis=1)
        pi1 = np.nansum(spi1) / nvalid[0] if nvalid[0] else np.nan
        pi2 = np.nansum(spi2) / nvalid[1] if nvalid[1] else np.nan
        pitot = np.nansum(spitot) / nvalid_tot if nvalid_tot else np.nan
        dxy = np.nansum(sdxy) / nvalid_tot if nvalid_tot else np.nan
        fst = (pitot - (pi1 + pi2)/2) / pitot if pitot else np.nan
        return pi1, pi2, pitot, dxy, fst, nvalid_tot


def read_mafs(
        mafs_files: list[Path],
        region: Region,
        pops: list[str],
        minind: int=0,
        minfreq: float=0.00001,
        interval: int=100000
        ) -> np.array:
    freqs = Frequencies(region, pops)
    for pidx, mafs_file in enumerate(mafs_files):
        logger.info(f"Working on file {mafs_file} ...")
        with TabixFile(str(mafs_file)) as mafs_in:
            try:
                header = mafs_in.header[0].split('\t')
            except Exception as e:
                logger.error(f"Couldn't open header of mafs file {mafs_file}!")
            processed = 0
            try: # Pysam throws a ValueError if the contig is not found in the tabix file
                for row in mafs_in.fetch(region.chrom, region.start-1, region.end, parser=asTuple()):
                    idx = int(row[1]) - region.start
                    nind = int(row[7])
                    freqs.nind[idx, pidx] = nind
                    if nind >= minind:
                        freq = DTYPE(row[6])
                        freq = 0. if freq < minfreq else freq
                        major, minor = Base[row[2]], Base[row[3]]
                        base_freqs = np.zeros(4)
                        base_freqs[major] = 1 - freq
                        base_freqs[minor] = freq
                        freqs.freqs[idx, pidx,:] = base_freqs
                    processed += 1
                    if not processed % interval: logger.info(f"Processed {processed} lines of mafs file.")
            except ValueError as e:
                logger.error(f"Couldn't find contig {region.chrom} in mafs file!")
                logger.error(e)
            logger.info(f"Processed {processed} lines of mafs file.")
    return freqs


class Frequencies_1b():
    def __init__(self, region: Region, pops: list[str]):
        self.freqs = np.empty((region.length, len(pops)), dtype=DTYPE)
        self.freqs.fill(np.nan)
        self.bases = np.zeros((region.length, 3), dtype=np.byte)
        self.nind = np.zeros((region.length, len(pops)), dtype=np.uint8)
        self.region = region
        self.pops = pops

    def is_valid(self, pos: int, pidx: int) -> bool:
        return np.isnan(self.freqs[pos-self.region.start, pidx])

    def get_position(self, pos: int, pidx: int) -> float:
        return self.freqs[pos-self.region.start, pidx]
    
    def get_position_row(self, pos: int) -> np.array:
        return self.freqs[pos-self.region.start,:]

    def set_position(self, pos: int, pidx: int, freq: float):
        self.freqs[pos-self.region.start, pidx] = freq

    def get_pop(self, spos: int, epos: int, pidx: int) -> np.array:
        sidx = spos-self.region.start
        eidx = epos-self.region.start + 1
        return self.freqs[sidx:eidx,pidx]
    
    def get_diversity(self, spos: int, epos: int) -> tuple[float,float,float,int]:
        sidx = spos-self.region.start
        eidx = epos-self.region.start + 1
        window = self.freqs[sidx:eidx,:]
        nvalid = np.sum(~np.isnan(window), 0)
        nvalid_tot = np.sum(~(np.isnan(window[:,0]) | np.isnan(window[:,1])))
        logger.info(f"nvalid pop1:{nvalid[0]}, nvalid pop2: {nvalid[1]}, nvalid total: {nvalid_tot}")
        spi1 = 2 * window[:,0] * (1 - window[:,0])
        spi2 = 2 * window[:,1] * (1 - window[:,1])
        sdxy = window[:,0] * (1 - window[:,1]) + window[:,1] * (1 - window[:,0])
        pi1 = np.nansum(spi1) / nvalid[0] if nvalid[0] else np.nan
        pi2 = np.nansum(spi2) / nvalid[1] if nvalid[1] else np.nan
        dxy = np.nansum(sdxy) / nvalid_tot if nvalid_tot else np.nan
        return pi1, pi2, dxy, nvalid_tot


def read_mafs_1b(
        mafs_files: list[Path],
        region: Region,
        pops: list[str],
        minind: int=0,
        minfreq: float=0.00001,
        interval: int=100000
        ) -> np.array:
    freqs = Frequencies_1b(region, pops)
    for pidx, mafs_file in enumerate(mafs_files):
        logger.info(f"Working on file {mafs_file} ...")
        with TabixFile(str(mafs_file)) as mafs_in:
            try:
                header = mafs_in.header[0].split('\t')
            except Exception as e:
                logger.error(f"Couldn't open header of mafs file {mafs_file}!")
            processed = 0
            try: # Pysam throws a ValueError if the contig is not found in the tabix file
                for row in mafs_in.fetch(region.chrom, region.start-1, region.end, parser=asTuple()):
                    idx = int(row[1]) - region.start
                    nind = int(row[7])
                    freqs.nind[idx, pidx] = nind
                    if nind >= minind:
                        freq = DTYPE(row[6])
                        freq = 0. if freq < minfreq else freq
                        bases = np.array([Base[base] for base in row[2:5]])
                        if freqs.bases[idx,0]: # major base is not N
                            if bases[0] != freqs.bases[idx,0]:
                                logger.critical(f"Major base doesn't match for position {int(row[1])}: {bases[0]} vs. {freqs.bases[idx,0]}!")
                        else:
                            freqs.bases[idx,0] = bases[0]
                        if freq > 0: # minor base only matters if freq is > 0
                            if freqs.bases[idx,1]: # minor base is not N
                                if bases[1] != freqs.bases[idx,1]:    # minor base only matters if freq is > 0
                                    logger.critical(f"Minor base doesn't match for position {int(row[1])}: {bases[1]} vs. {freqs.bases[idx,1]}!")
                            else:
                                freqs.bases[idx,1] = bases[1]
                        freqs.freqs[idx, pidx] = freq
                    processed += 1
                    if not processed % interval: logger.info(f"Processed {processed} lines of mafs file.")
            except ValueError as e:
                logger.error(f"Couldn't find contig {region.chrom} in mafs file!")
                logger.error(e)
            logger.info(f"Processed {processed} lines of mafs file.")
    return freqs


def run_test():
    logging.basicConfig(level="INFO", format="[%(levelname)s] %(message)s")

    pops = ["Izinga_blue", "Izinga_yellow"]
    mafs = [Path("/data/users/anater/pipelines/nf-core-popgen/results_test_le5mb/angsd/Izinga_blue.mafs.gz"), \
            Path("/data/users/anater/pipelines/nf-core-popgen/results_test_le5mb/angsd/Izinga_yellow.mafs.gz")]
    region = Region.from_string("tarseq_0_arrow:1000000-2000000")

    freqs = read_mafs_4b(mafs, region, pops)
    div = freqs.get_diversity(1200000, 1400000)
    print(div)

    freqs2 = read_mafs(mafs, region, pops)
    div2 = freqs2.get_diversity(1200000, 1400000)
    print(div2)
