#!/usr/bin/python
# -*- coding: utf-8 -*-
import eospac as eos
import numpy as np
from numpy.testing import assert_allclose

from scipy.constants import physical_constants
R_CST = physical_constants['molar gas constant'][0]*1e7 # erg.K⁻¹.mol⁻¹



def setup():
    global tables_list, material
    global eosmat
    global abar, gamc, game


    tables_list = ['.*_DT', '.*_DUt']
    material = None
    gamc, game, abar = 1.66, 1.66, 3.5
    eosmat = eos.EosMaterial(material,
            tables=tables_list,
            options={'game': game, 'gamc': gamc, 'abar': abar},
            backend='gamma')

def test_interpolation():
    """Checking the gamma law partial derivatives"""
    #for table in tables_list:
    #    ctable = getattr(sestab, table)
    #    R_array = ctable['R_Array']
    #    T_array = ctable['T_Array']
    #    F_array = ctable['F_Array']

    #    R,T = np.meshgrid(R_array, T_array, indexing='ij')

    #    F_array_itp = ctable(R, T)
    #    #F_array_itp = F_array_itp.T.reshape(F_array.shape)
    #    yield np.testoing.assert_allclose, F_array, F_array_itp, 1e-15
    Npt  = 10
    rho = np.logspace(-5, 1, Npt)
    temp = np.linspace(100, 10000, Npt)

    yield assert_allclose, eosmat.Pt_DT(rho, temp), rho*R_CST*temp/abar
    yield assert_allclose, eosmat.Pt_DT.dFx(rho, temp), R_CST*temp/abar
    yield assert_allclose,   eosmat.Pt_DT.dFxx(rho, temp), np.zeros(Npt), 1e-8, 20
    yield assert_allclose,   eosmat.Pt_DT.dFxy(rho, temp), R_CST/abar*np.ones(Npt)
    yield assert_allclose, eosmat.Pt_DT.dFy(rho, temp), rho*R_CST/abar

#print al_eos.Ut_DT(X,Y)
#for tab in all_tables:
#    print tab, getattr(al_eos, tab).options
#print al_eos.Ut_DT


