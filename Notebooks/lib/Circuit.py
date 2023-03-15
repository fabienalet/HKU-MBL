#!/usr/bin/env python
from scipy.sparse import csr_matrix, kron
import numpy as np
from scipy.linalg import svd

"""
Code by Nicolas Macé (mace@irsamc.ups-tlse.fr)
"""

def diagu(U):
    """
    Diagonalize a unitary matrix, taking advantage
    of the fast diagonalization of the Hermitian matrix H = U + U.H
    """
    from scipy.linalg.lapack import zheevd
    H = U + U.H
    diagH, P, info = zheevd(H)
    if info:
        raise RuntimeError(f'zheevd failed with error code {info}')
    # check for degeneracies
    tol = 1e-14
    ndegen = len((np.abs(diagH[1:] - diagH[:-1]) < tol).nonzero()[0])
    assert ndegen == 0, 'Degeneracies in the spectrum of H = U + U.H.'
    # recover the eigenvalues of U
    diagU = np.diag(P.T.conj() @ U @ P)
    return diagU, P

def sparse_id(N):
    """
    Construct the identity as a sparse matrix
    """
    r = range(N)
    return csr_matrix(([1.]*N, (r, r)))

def sparse_diag(diag):
    r = range(len(diag))
    return csr_matrix((diag, (r, r)))

class KickedIsing():
    """
    Kicked Ising model.

    One recovers the model defined in https://arxiv.org/abs/1612.07327 (eqns 31-32) doing
    -- hx = [tau*hx]*L
    -- hy = [tau*hy]*L
    -- hz = [tau*hz]*L
    -- Jz = [tau]*L
    where the LHS is in the notation of the article.

    One recovers the model defined in https://arxiv.org/abs/1608.02589 doing
    -- hx = [g - epsilon]*L
    -- hy = [0]*L
    -- hz = [Bz[i] for i in range(L)]
    -- Jz = Jz
    where the LHS is in the notation of the article.

    """
    def __init__(self, L, hx, hy, hz, Jz):
        self.L = L
        self.d = 2 # we work with qubits
        self.hx = hx.copy()
        self.hy = hy.copy()
        self.hz = hz.copy()
        self.Jz = Jz.copy()
        self.compute_XY()
        self.compute_Z()
        
        self._Z = []
        for r in range(L):
            self._Z.append(self._local_Z_operator(r))

    def _local_XY_unitary(self, r):
        """
        Return the X + Y unitary in the space of qubit at r
        """
        nleft = r # number of qubits to the left of qubit at r
        nright = self.L-r-1 # number of qubits to the right of qubit at r+1
        X = np.array([[0.,1.],[1.,0.]])
        Y = np.array([[0.,-1j],[1j,0.]])
        XYterm = self.hx[r] * X + self.hy[r] * Y
        # P @ np.diag(diag) @ np.conjugate(P).T = XYterm
        diag, P = np.linalg.eig(XYterm)
        # exp(-i*(hx X + hy Y))
        mat = P @ np.diag(np.exp(-1j*diag)) @ np.conjugate(P).T
        rot = csr_matrix(mat)
        if nleft > 0:
            rot = kron(sparse_id(self.d**nleft), rot)
        if nright > 0:
            rot = kron(rot, sparse_id(self.d**nright))
        return rot

    def compute_XY(self):
        """
        Rotation by fields along the X and Y directions, stored as L sparse matrices
        """
        self.UXYs = []
        for r in range(self.L):
            self.UXYs.append(self._local_XY_unitary(r))

    def _local_ZZ_unitary(self, r):
        """
        Return the ZZ unitary in the space of qubits at r and r+1
        """
        nleft = r # number of qubits to the left of qubit at r
        nright = self.L-r-2 # number of qubits to the right of qubit at r+1

        # interaction term in the (00, 01, 10, 11) basis
        ZZ = self.Jz[r] * np.array([1., -1, -1, 1])
        ZZ = sparse_diag(np.exp(-1j*ZZ))
        if nleft > 0:
            ZZ = kron(sparse_id(self.d**nleft), ZZ)
        if nright > 0:
            ZZ = kron(ZZ, sparse_id(self.d**nright))
        return ZZ

    def _local_Z_unitary(self, r):
        """
        Return the Z unitary in the space of qubit at r
        """
        nleft = r # number of qubits to the left of qubit at r
        nright = self.L-r-1 # number of qubits to the right of qubit at r+1

        # field term
        Z = np.array([self.hz[r], -self.hz[r]])
        Z = sparse_diag(np.exp(-1j*Z))
        if nleft > 0:
            Z = kron(sparse_id(self.d**nleft), Z)
        if nright > 0:
            Z = kron(Z, sparse_id(self.d**nright))
        return Z

    def _local_Z_operator(self, r):
        """
        Return the Z Pauli operator at qubit r
        """
        nleft = r # number of qubits to the left of qubit at r
        nright = self.L-r-1 # number of qubits to the right of qubit at r+1

        # field term
        Z = np.array([1., -1.])
        Z = sparse_diag(Z)
        if nleft > 0:
            Z = kron(sparse_id(self.d**nleft), Z)
        if nright > 0:
            Z = kron(Z, sparse_id(self.d**nright))
        return Z

    def compute_Z(self):
        """
        Nearest neighbors dephasing + rotation along Z
        """
        self.nconf = 2**self.L
        self.UZ = sparse_id(self.nconf)
        # interaction part (open bc)
        for r in range(self.L-1):
            self.UZ *= self._local_ZZ_unitary(r)
        # field part
        for r in range(self.L):
            self.UZ *= self._local_Z_unitary(r)
            
    def Z(self, r):
        """
        Return the Z operator at site r
        """
        return self._Z[r]
            
    def apply_Z(self, r, vec):
        """
        Return sigmaz(r)*vec
        TODO
        """
        return self._Z[r] @ vec

    def apply(self, vec):
        """
        matrix-vector product doesn't use OMP threads. Use SpMV_viaMKL for that.
        """
        vec = self.UZ.dot(vec)
        for UXY in self.UXYs:
            vec = UXY.dot(vec)
        return vec

    def todense(self):
        """
        Convert to dense (numpy) matrix
        """
        U = self.UZ.todense() # way faster than performing sparse mat-mat product
        for UXY in self.UXYs:
            U = UXY @ U
        return U

if __name__ == "__main__":
    L = 20
    nconf = 2**L
    LA = L//2
    tau = 0.8
    hx = [tau*0.9045]*L
    hy = [tau*0.3457]*L
    hz = [tau*0.8090]*L
    Jz = [tau]*L
    op = KickedIsing(L, hx, hy, hz, Jz)
    psi = np.zeros(nconf)
    i0 = np.random.randint(0, nconf)
    psi[i0] = 1.
    nconfB = 2**(L-LA)
    for t in range(10):
        spec = entanglement_spectrum(nconfB, psi)
        ent = np.sum(H(spec, 1e-16))
        psi = op.apply(psi)
        print(f'S({t})={ent}')
