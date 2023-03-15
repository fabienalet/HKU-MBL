"""
Code by Nicolas Macé (mace@irsamc.ups-tlse.fr)
"""

class Ising_conf():
    def __init__(self, qubits: str):
        self._qubits = qubits
        
    def __add__(self, conf):
        return Ising_conf(self._qubits + conf._qubits)
    
    def __str__(self):
        return self._qubits
    
    def __repr__(self):
        return self._qubits
    
    def __hash__(self):
        return int(self._qubits, 2)
    
class Ising_basis_generator():
    """
    Store a configuration and generate the next one
    """
    def __init__(self, L: int):
        self._L = L
        self._format = f'0{self._L}b'
        self._hash = 0
        self._conf = Ising_conf(format(self._hash, self._format))
        self._maxhash = 2**self._L
        
    def __next__(self):
        if self._hash == self._maxhash:
            raise StopIteration
        else:
            self._hash += 1
            old_conf = self._conf
            self._conf = Ising_conf(format(self._hash, self._format))
            return old_conf
            
    def __iter__(self):
        return self
    
class IsingBasis:
    """
    Mimicks Basis.h
    """
    def __init__(self, L):
        # number of qubits
        self._L = L
        # 
        self._format = f'0{self._L}b'
        # list of confs (cast to tuple to have immutable unique objects)
        self._confs = [conf for conf in Ising_basis_generator(self._L)]
        self._InverseMap = {hash(conf):i for i, conf in enumerate(self._confs)}

    def index(self, conf):
        """
        Return index of the conf
        """
        return self._InverseMap[hash(conf)]
    
class IsingOperator():
    def __init__(self, arg, basis: IsingBasis):
        (data, (row_conf, col_conf)) = arg
        self._L = basis._L
        self._basis = basis
        self._row_conf = row_conf
        self._col_conf = col_conf
        self._data = data
        row = [basis.index(cur_row_conf) for cur_row_conf in row_conf]
        col = [basis.index(cur_col_conf) for cur_col_conf in col_conf]
        
    def to_coo(self):
        row = [self._basis.index(row_conf) for row_conf in self._row_conf]
        col = [self._basis.index(col_conf) for col_conf in self._col_conf]
        size = 2**self._L
        return sparse.coo_matrix((self._data, (row, col)), shape=(size, size))
    
def tensor(a: IsingOperator, b: IsingOperator):
    """
    Tensor product of a and b
    """
    basis = IsingBasis(a._L + b._L)
    data, row, col = [], [], []
    for data_a, row_conf_a, col_conf_a in zip(a._data, a._row_conf, a._col_conf):
        for data_b, row_conf_b, col_conf_b in zip(b._data, b._row_conf, b._col_conf):
            data.append(data_a * data_b)
            row.append(row_conf_a + row_conf_b)
            col.append(col_conf_a + col_conf_b)
    return IsingOperator((data, (row, col)), basis)