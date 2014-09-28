#!/usr/bin/python # -*- coding: utf-8 -*-
import eospac as eos
import numpy as np
from numpy.testing import assert_allclose

from scipy.constants import physical_constants
R_CST = physical_constants['molar gas constant'][0]*1e7 # erg.K⁻¹.mol⁻¹



def setup():
    global tables_list, material
    global eosmat
    global abar, gamc, game


    tables_list = ['St_DUt', 'Pt_DT', 'T_DUt', 'Ut_DT', 'At_DT', 'Pt_DUt', 'T_DPt', 'Pt_DSt', 'St_DT']
    # http://courses.chem.indiana.edu/c360/documents/vanderwaalsconstants.pdf
    material = None
    # parameters for water
    a =  5.537 # bar*L²/mol² 
    b = 0.0305 # L/mol
    delta = 0.329
    abar = 6.0
    eosmat = eos.EosMaterial(material,
            tables=tables_list,
            options={'delta': delta, 'a': a, 'b': b, 'abar': abar},
            units='eospac',
            backend='vdw')

def test_critical_point():
    rho_c = 0.0655737704918 # g/cm³
    Pt_c = 22.06e-3 # GPa
    temp_c = 647.096 # K

    yield assert_allclose, eosmat.Pt_DT['Pt_c']*1e-10, Pt_c, 1e-2
    yield assert_allclose, eosmat.Pt_DT['rho_c'], rho_c, 1e-2
    yield assert_allclose, eosmat.Pt_DT['temp_c'], temp_c, 1e-2
    yield assert_allclose, eosmat.Pt_DT['Pt_c']/(eosmat.Pt_DT['temp_c']*eosmat.Pt_DT['rho_c']), 3*R_CST/(8*eosmat.Pt_DT['abar'])
    yield assert_allclose, Pt_c*1e10/(rho_c*temp_c), 3*R_CST/(8*eosmat.Pt_DT['abar']), 1e-3
    yield assert_allclose, eosmat.Pt_DT(rho_c, temp_c), eosmat.Pt_DT['Pt_c']*1e-10, 1e-3

def test_inverse_tabs():
    rho = 0.05
    temp = 1000
    pres = eosmat.Pt_DT(rho, temp)
    eint = eosmat.Ut_DT(rho, temp)
    sint = eosmat.St_DT(rho, temp)
    aint = eosmat.At_DT(rho, temp)
    yield assert_allclose,  eosmat.St_DUt(rho, eint), sint, 1e-07, 0, "St_DUt"
    yield assert_allclose,  eosmat.T_DUt(rho, eint), temp, 1e-7, 0, 'T_DUt'
    yield assert_allclose,  eosmat.Pt_DUt(rho, eint), pres, 1e-7, 0, 'Pt_DUt'
    yield assert_allclose,  eosmat.T_DPt(rho, pres), temp, 1e-7, 0, 'T_DPt'



def test_fund_derivative():
    # Comparing with values in Guardone et al 2001
    rho_c = eosmat.Pt_DT['rho_c']
    Pt_c = eosmat.Pt_DT['Pt_c']
    rho = np.array([1.01, 0.594])*rho_c
    pres = np.array([1.6077, 0.89570])*Pt_c
    temp = eosmat.T_DPt(rho, pres)
    G3 = np.array([2.23946, 1.36136])
    #for idx in range(2):
        #yield  assert_allclose, eosmat.q['G3_vdw'](rho[idx], temp[idx]), G3[idx], 1e-4

def test_consitency_wdv_ideal_gas():
    rho = 1e-6
    temp = 10*11640
    eosmat_gamma = eos.EosMaterial(material,
              tables=['Pt_DT'],
              options={'game': 1.6666, 'abar': abar, 'gamc': 1.6666},
              units='eospac',
              backend='gamma')
    yield  assert_allclose, eosmat.Pt_DT(rho, temp), eosmat_gamma.Pt_DT(rho, temp), 1e-5



def test_entropy_validity():
    rho_c = eosmat.Pt_DT['rho_c']
    Pt_c = eosmat.Pt_DT['Pt_c']*1e-10
    temp_c = eosmat.Pt_DT['temp_c']
    eint_c = eosmat.Ut_DT(rho_c, temp_c)
    sint_c = eosmat.St_DUt(rho_c, eint_c)
    yield assert_allclose, eosmat.T_DUt(rho_c, eint_c), temp_c
    yield assert_allclose, eosmat.Pt_DSt(rho_c, sint_c), Pt_c


