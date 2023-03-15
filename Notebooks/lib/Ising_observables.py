#!/usr/bin/env python
import numpy as np
from scipy.linalg import svd

"""
Code by Nicolas Macé (mace@irsamc.ups-tlse.fr)
"""

"""
Observables in the case of a model with spin-1/2 variables (eg Ising)
"""

H = lambda x, cutoff: np.where(x < cutoff, 0., -x*np.log(x))

def entanglement_spectrum(LA: int, state: np.array):
    """
    Entanglement spectrum tracing over A
    """
    nconfA = 2**LA
    # column index varies on one subsystem, row on the other
    state = state.reshape(nconfA, -1).copy()
    # singular values are the square roots of the entanglement eigenvalues
    return svd(state, compute_uv=False, overwrite_a=True)**2.

def vNM_entropy(LA: int, state: np.array):
    """
    von Neumann entropy
    """
    cutoff = 1e-16
    es = entanglement_spectrum(LA, state)
    return np.sum(H(es, cutoff))

def rgaps(spec):
    """
    Compute the gap ratios
    """
    # sort
    spec = np.sort(np.nan_to_num(spec))
    # level spacings
    spac = (np.roll(spec,-1) - spec)[:-1]
    # rgaps
    rgap = (np.roll(spac,-1)/spac)[:-1]
    # normalized rgaps
    rgap = np.min([rgap, 1./rgap], axis=0)
    return rgap
