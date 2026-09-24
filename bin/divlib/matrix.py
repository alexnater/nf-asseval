import logging
import array as arr

logger = logging.getLogger(__name__)

class Matrix:
    def __init__(self, type: str, nrows: int=0, ncols: int=0, init=0):
        self.data = arr.array(type, [init] * nrows * ncols)
        self.nrows = nrows
        self.ncols = ncols
        self.max = (1 << (self.data.itemsize * 8)) - 1

    def __getitem__(self, pos: "tuple[int,int]"):
        return self.data[pos[0] * self.ncols + pos[1]]

    @property
    def shape(self):
        return self.nrows, self.ncols
    
    @property
    def length(self):
        return len(self.data)
    
    def add_rows(self, nrows: int, init=0) -> int:
        self.data.extend([init] * nrows * self.ncols)
        self.nrows += nrows
        return self.nrows

    def add_set_row(self, values: list):
        if len(values) != self.ncols:
            raise IndexError("length of values doesn't correspond to number of columns in matrix")
        self.data.extend(values)
        self.nrows += 1

    def get_row(self, row: int) -> arr.array:
        if row >= self.nrows:
            raise IndexError("array index out of range")
        start = row * self.ncols
        return self.data[start:(start + self.ncols)]

    def set_row(self, row: int, values):
        if row >= self.nrows:
            raise IndexError("array index out of range")
        start = row * self.ncols
        self.data[start:(start + self.ncols)] = arr.array(self.data.typecode, values)

    def get_col(self, col: int) -> arr.array:
        if col >= self.ncols:
            raise IndexError("array index out of range")
        return arr.array(self.data.typecode, self.data[col::self.ncols])
    
    def set_col(self, col: int, values):
        if col >= self.ncols:
            raise IndexError("array index out of range")
        self.data[col::self.ncols] = arr.array(self.data.typecode, values)

    def get(self, row: int, col: int):
        if row >= self.nrows or col >= self.ncols:
            raise IndexError("array index out of range")
        return self.data[row * self.ncols + col]
    
    def set(self, row: int, col: int, value):
        if row >= self.nrows or col >= self.ncols:
            raise IndexError("array index out of range")
        self.data[row * self.ncols + col] = value
