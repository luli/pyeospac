#!/usr/bin/python
# -*- coding: utf-8 -*-

from eospac.base import EosUnits
from eospac import EosMaterial
from eospac.eospac.constants import tables

import numpy as np
from eospac.base import eV2K_cst

from numpy.testing import assert_allclose

def is_equal(x,y):
    assert x == y

def test_consistency():
    """ Check consistency"""

    u1 = EosUnits('cgs', 'cgs')
    yield is_equal, u1.o2r('P'), 1.

def test_table_parsing():
    """ Checking that the tables are parsed correctly """
    u1 = EosUnits('eospac', 'cgs')
    for tab in tables:
        # exculding fake tables and non implemented EOSPAC tables
        if '_' in tab and np.all([sym not in tab for sym in ['K', 'G', 'B', 'M']]):
            u1.table_units(tab)

def test_units_gamma():
    """ Test units conversion with gamma backend """
    unit_systems = ['cgs', 'eospac', 'feos']
    T_conv = dict(zip(unit_systems,(1.0, 1.0, 1./eV2K_cst)))
    U_conv = dict(zip(unit_systems,(1.0, 1e10, 1)))

    meint, meint_dx, meint_dy, meint_dxy = {}, {}, {}, {}
    
    for unit_system in unit_systems:
        mat = EosMaterial(backend='gamma', tables=['.*_D.*'],
                            units=unit_system,
                            options={'game':1.66, 'abar': 3.5, 'gamc': 1.66})
        rho = 1.04
        tele = 300*T_conv[unit_system]
        meint[unit_system] = mat.Ut_DT(rho, tele)
        meint_dx[unit_system] = mat.Ut_DT.dFx(rho, tele)
        meint_dy[unit_system] = mat.Ut_DT.dFy(rho, tele)
        meint_dxy[unit_system] = mat.Ut_DT.dFxy(rho, tele)
    us0 = unit_systems[0]
    for usi in unit_systems[1:]:
        yield assert_allclose, meint[us0], meint[usi]*U_conv[usi]
        yield assert_allclose, meint_dx[us0], meint_dx[usi]*U_conv[usi]
        yield assert_allclose, meint_dy[us0], meint_dy[usi]*U_conv[usi]*T_conv[usi]
        yield assert_allclose, meint_dxy[us0], meint_dxy[usi]*U_conv[usi]*T_conv[usi]

def test_compare_eospac_gamma():
    """Check that sesame and gamma give simular results in the ideal gas region"""
    rho = 1e-3
    tele = 1000.*eV2K_cst
    mat_g = EosMaterial(backend='gamma', tables=['Pt_DT'],
                            units='eospac',
                            options={'game':1.26, 'abar': 3.5, 'gamc': 1.66})
    mat_s = EosMaterial(backend='eospac', tables=['Pt_DT'],
                            units='eospac',
                            material=7593,
                            options={})
    #yield assert_allclose, mat_g.Pt_DT(rho, tele), mat_s.Pt_DT(rho, tele)
    #yield assert_allclose, mat_g.q['Ct2'](rho, tele)**0.5, mat_s.q['Ct2'](rho, tele)**0.5












