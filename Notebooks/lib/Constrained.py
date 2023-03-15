import scipy.sparse as sparse
import numpy as np
from Constrained_basis import Conf, Basis

"""
Code by Nicolas Macé (mace@irsamc.ups-tlse.fr)
"""


class Operator():
    def __init__(self, arg, basis):
        (data, (row_conf, col_conf)) = arg
        self._L = basis._L
        self._basis = basis
        self._row_conf = row_conf
        self._col_conf = col_conf
        self._data = data

    @property
    def L(self):
        return self._L

    @property
    def size(self):
        return len(self._data)

    @property
    def data(self):
        return self._data

    @property
    def row_conf(self):
        return self._row_conf

    @property
    def col_conf(self):
        return self._col_conf

    def tocoo(self):
        row = [self._basis.index(row_conf) for row_conf in self._row_conf]
        col = [self._basis.index(col_conf) for col_conf in self._col_conf]
        size = self._basis.Nobc
        return sparse.coo_matrix((self._data, (row, col)), shape=(size, size))

    def __repr__(self):
        return self.tocoo().todense().__repr__()

def kron(a: Operator, b: Operator):
    """
    Tensor product of a and b
    """
    as_, bs = a.size, b.size
    if as_ == 0 or bs == 0:
        if bs != 0:
            return b
        elif as_ != 0:
            return a
        else:
            raise RuntimeError('Attempting Kronecker product of empty operators')
    basis = Basis(a.L + b.L)
    data, row, col = [], [], []
    for data_a, row_conf_a, col_conf_a in zip(a.data, a.row_conf, a.col_conf):
        for data_b, row_conf_b, col_conf_b in zip(b.data, b.row_conf, b.col_conf):
            try:
                cur_row_conf = row_conf_a + row_conf_b
                cur_col_conf = col_conf_a + col_conf_b
                row.append(cur_row_conf)
                col.append(cur_col_conf)
                data.append(data_a * data_b)
            except KeyError:
                pass
    return Operator((data, (row, col)), basis)
