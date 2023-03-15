#!/usr/bin/env python

"""
Code by Nicolas Macé (mace@irsamc.ups-tlse.fr)
"""

def iterative_OTOC(H, psi_0, times):
    """
    Iterative version of the OTOC evaluation: we do not compute the evolution operator explicitly, but using Krylov techniques.
    This is faster than full diagonalization for large enough systems and small enough times.
    """
    # evolution operator
    U = exp_op(H,a=-1j,start=times.min(),stop=times.max(),num=len(times),iterate=True)
    OTOC = []
    cur_steps = 1
    for cur_t, Upsi, UVpsi in tqdm(zip(times, U.dot(psi_0), U.dot(V.dot(psi_0)))):
        ## OTOC
        WbUpsi, WUVpsi = W.H.dot(Upsi), W.dot(UVpsi)
        # backward evolution operator
        Ub = exp_op(H,a=1j,start=times.min(),stop=cur_t,num=cur_steps,iterate=True)
        for UbWbUpsi, UbWUVpsi in zip(Ub.dot(WbUpsi), Ub.dot(WUVpsi)):
            pass
        VUbWUVpsi = V.dot(UbWUVpsi)
        cur_OTOC = UbWbUpsi.conj() @ VUbWUVpsi
        OTOC.append(cur_OTOC)
        cur_steps += 1
    return OTOC