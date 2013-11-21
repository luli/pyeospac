#!/usr/bin/python
# -*- coding: utf-8 -*-
import eospac as eos
import numpy as np
import os.path
from numpy.testing import assert_allclose

from scipy.constants import physical_constants
R_CST = physical_constants['molar gas constant'][0]*1e7 # erg.K⁻¹.mol⁻¹



def test_ionmix():
    """Reading an ionmix table"""
    mpath = "/home/rth/luli/NestedOutflows/NestedOutflows/simulations/NestedOutflows"
    mat1 = eos.EosMaterial(63720,
            options={'type': 'ionmix',
                'abar':  26.9815, 'zbar':13.0, 'rho_ref': 2.7,
                    'path': os.path.join(mpath, 'al-imx-32g.cn4')},
            units='cgs', backend='tabulated')
    mat1.save('/tmp/ionmix_63720.ses')

    mat2 = eos.EosMaterial(63720, tables=['.t_DT'],
            units='cgs', backend='eospac')

    temp = np.array([2.000E+03*eos.eV2K_cst])
    rho = np.array([4.480E-01])
    # comparing with values in al-imx-32g.imx
    yield assert_allclose, mat2.Pt_DT(rho, temp), np.array([3.204E+13+4.102E+14]), 1e-3
    yield assert_allclose, mat2.Ut_DT(rho, temp), np.array([1.188E+07+1.591E+08])*1e7, 1e-4



#    yield assert_allclose, eosmat.Pt_DT(rho, temp), rho*R_CST*temp/abar
#    yield assert_allclose, eosmat.Pt_DT.dFx(rho, temp), R_CST*temp/abar
#    yield assert_allclose,   eosmat.Pt_DT.dFxx(rho, temp), np.zeros(Npt), 1e-8, 20
#    yield assert_allclose,   eosmat.Pt_DT.dFxy(rho, temp), R_CST/abar*np.ones(Npt)
#    yield assert_allclose, eosmat.Pt_DT.dFy(rho, temp), rho*R_CST/abar

#print al_eos.Ut_DT(X,Y)
#for tab in all_tables:
#    print tab, getattr(al_eos, tab).options
#print al_eos.Ut_DT


