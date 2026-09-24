import sys
import logging
from pathlib import Path
import numpy as np
from pysam import TabixFile, asTuple
from divlib.general import Base
from divlib.region import Region

logger = logging.getLogger(__name__)

DTYPE = np.uint16
MAX = np.iinfo(DTYPE).max


class Depth():
    def __init__(self, region: Region, nsamples: int, pops: np.array=None):
        self.depth = np.zeros((region.length, nsamples), dtype=DTYPE)
        self.region = region
        self.pops = pops

    def get_indices(self, start_pos: int, end_pos: int=None) -> tuple[int,int]:
        start_idx = start_pos - self.region.start
        end_idx = end_pos - self.region.start + 1 if end_pos else start_idx + 1
        return start_idx, end_idx

    def is_valid(self, pos: int, sidx: int, mindepth: int) -> bool:
        return self.depth[pos-self.region.start, sidx] >= mindepth

    def is_valid_row(self, pos: int, mindepth: int) -> np.array:
        return np.array(self.depth[pos-self.region.start,:] >= mindepth, dtype=np.bool)

    def get_position(self, pos: int, sidx: int):
        return self.depth[pos-self.region.start, sidx]
    
    def get_position_row(self, pos: int) -> np.array:
        return np.array(self.depth[pos-self.region.start,:], dtype=DTYPE)

    def set_position(self, pos: int, sidx: int, depth: int):
        self.depth[pos-self.region.start, sidx] = depth if depth <= MAX else MAX

    def set_row(self, pos: int, depths: list[int]):
        row = np.array(tuple(dp if dp <= MAX else MAX for dp in depths), dtype=DTYPE)
        self.depth[pos-self.region.start,:] = row

    def get_colsum(self, start_pos: int, end_pos: int, sidx: int) -> np.uint64:
        start_idx, end_idx = self.get_indices(start_pos, end_pos)
        return np.sum(self.depth[start_idx:end_idx, sidx], dtype=np.uint64)

    def get_colsum_row(self, start_pos: int, end_pos: int) -> np.array:
        start_idx, end_idx = self.get_indices(start_pos, end_pos)
        return np.sum(self.depth[start_idx:end_idx,:], axis=0, dtype=np.uint64)

    def get_nvalid(self, start_pos: int, end_pos: int, sidx: int, mindepth: int) -> np.uint32:
        start_idx, end_idx = self.get_indices(start_pos, end_pos)
        return np.sum(self.depth[start_idx:end_idx, sidx] >= mindepth, dtype=np.uint32)

    def get_nvalid_row(self, start_pos: int, end_pos: int, mindepth: int) -> np.array:
        start_idx, end_idx = self.get_indices(start_pos, end_pos)
        return np.sum(self.depth[start_idx:end_idx,:] >= mindepth, axis=0, dtype=np.uint32)


def read_depth(
        depth_file: Path,
        region: Region,
        samples: "list[str]",
        interval: int=1000
        ) -> Depth:
    if not region.defined:
        raise ValueError("read_depth_array requires a valid region!")
    depths = Depth(region, len(samples))
    with TabixFile(str(depth_file)) as depth_in:
        try:
            header = depth_in.header[0].split('\t')
            columns_by_id = get_column_dict(header, samples)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        row_idx = np.array([columns_by_id[sample] for sample in samples], dtype=np.uint16)
        processed = 0
        try: # Pysam throws a ValueError if the contig is not found in the tabix file
            for row in depth_in.fetch(region.chrom, region.start-1, region.end, parser=asTuple()):
                depths.depth[int(row[1])-region.start,:] = np.array(row)[row_idx]
                processed += 1
                if not processed % interval: logger.info(f"Processed {processed} lines of depth file.")
        except ValueError as e:
            logger.error(f"Couldn't find contig {region.chrom} in depth file!")
        logger.info(f"Processed {processed} lines of depth file.")
    return depths


def get_column_dict(header: "list[str]", samples: "list[str]", validate: bool=True) -> dict:
    columns_by_id = {Path(fn).with_suffix('').name.replace('_A', ''): col for col, fn in enumerate(header) if col >= 2}
    for sample in samples:
        if sample in columns_by_id:
            logger.info(f"Found {sample} at column {columns_by_id[sample] + 1} in depth file.")
        elif validate:
            raise Exception(f"Sample {sample} not found in depth file!")
        else:
            logger.warning(f"Sample {sample} not found in depth file!")
    return columns_by_id

def get_samples_from_depth_file(depth_file: Path) -> "list[str]":
    with TabixFile(str(depth_file)) as depth_in:
        header = depth_in.header[0].split('\t')
    return [Path(fn).name.split('.')[0] for fn in header[2:]]
