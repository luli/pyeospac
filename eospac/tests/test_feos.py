#!/usr/bin/python
# -*- coding: utf-8 -*-

from eospac.base import EosUnits
from eospac import EosMaterial
from eospac.feos import FeosMaterial
import eospac.feos.interface

import numpy as np
from eospac.base import eV2K_cst

from numpy.testing import assert_allclose

def is_equal(x,y):
    assert x == y

def test_low_level():
    """ Check consistency"""
    maxwell_temps = np.logspace(-4, 1, 50) 
    rho = np.logspace(-4, 1, 10)
    temp = np.logspace(-3, 3, 10)

    use_maxwell = 1
    table_handle = 1
    use_softsphere = 1
    matid = 3717
    my_opts ={'use_maxwell': True, 'use_softspheres': True,
                'maxwell_temp_arr': None, 'max_materials': 2,
                'grid_subsample_default': 0}
    mat1 = FeosMaterial(3717, tables=['.*_DT'], options=my_opts, units='cgs')
    D_arr, T_arr = mat1.Pt_DT['D_Array'], mat1.Pt_DT['T_Array']
    D, T = np.meshgrid(D_arr, T_arr, indices='ij')
    #print mat1.Pt_DT(D, T).min()
    #print mat1.Ut_DT(D, T).min()

    #Nel, N_maxwell_iso, Maxwell_T_reliable = _init_mat(table_handle, matid, use_maxwell, use_softsphere, maxwell_temps)
    #res = _get_eos_all(table_handle, Nel, use_maxwell, rho, temp)
    #res2 = _get_eos(table_handle, 0,  Nel, use_maxwell, rho, temp)
    #yield assert_allclose, res['Pt'], res2['P']
    #print _get_mat_par(table_handle, Nel)
    #print _get_crit_point(table_handle)
    #print _get_soft_shere_par(table_handle)

