from functools import lru_cache
import numpy as np

"""
Code by Nicolas Macé (mace@irsamc.ups-tlse.fr)
"""

class Conf():
    """
    Stores a valid configuration of the constrained Hilbert space
    """
    def __init__(self, qubits: str):
        if '11' in qubits:
            raise KeyError(f'Invalid configuration {qubits}')
        self._qubits = qubits

    def __add__(self, conf):
        return Conf(self._qubits + conf._qubits)

    def __str__(self):
        return self._qubits

    def __repr__(self):
        return self._qubits

    def __hash__(self):
        return int(self._qubits, 2)

zero = Conf('0')
one = Conf('1')
def next_sectors(S00, S01, S10, S11):
    """
    Return sectors obtained by adding one spin to the right of the input sectors
    """
    return [C + zero for C in S00 + S01], [C + one for C in S00], [C + zero for C in S10 + S11], [C + one for C in S10]

@lru_cache(maxsize=128)
def sectors(L: int):
    """
    Return sectors of the system of L spins
    """
    if L < 2:
        raise RuntimeError(f'Sectors are ill-defined for a chain of size L = {L}.')
    S00 = [zero + zero]
    S01 = [zero + one]
    S10 = [one + zero]
    S11 = []
    for _ in range(L-2):
        S00, S01, S10, S11 = next_sectors(S00, S01, S10, S11)
    return S00, S01, S10, S11

def basis_vectors(L: int):
    """
    Return a list of basis vectors of constrained system of L spins, using open boundary conditions
    """
    if L == 0:
        return []
    if L == 1:
        return [zero, one]
    return np.concatenate(sectors(L))

def size11(L: int):
    """
    Return the size of the 11 sector
    """
    if L == 0 or L == 1:
        return 0
    return len(sectors(L)[3])

class Basis:
    """
    List all configurations on the constrained chain of L sites
    """
    def __init__(self, L: int):
        # number of qubits
        self._L = L
        # list of confs
        self._confs = basis_vectors(L)
        # Hilbert space size for OBC
        self._Nobc = len(self._confs)
        # Hilbert space size for PBC
        self._Npbc = self._Nobc - size11(L)
        # index of the configurations
        self._InverseMap = {hash(conf):i for i, conf in enumerate(self._confs)}

    @property
    def L(self):
        return self._L

    @property
    def Nobc(self):
        return self._Nobc

    @property
    def Npbc(self):
        return self._Npbc

    @property
    def confs(self):
        return self._confs

    def index(self, conf):
        """
        Return index of the conf
        """
        return self._InverseMap[hash(conf)]
