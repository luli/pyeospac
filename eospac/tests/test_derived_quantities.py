#!/usr/bin/python
# -*- coding: utf-8 -*-
import eospac as eos
import numpy as np
import nose


def setup():
    global tables_list, material
    global tab
    global X, Y

    # Load the tables to compute all the derived quantities for all species
    material = 7593
    tab = eos.EosMaterial(material,
            eos.derived_quantities.keys(),
            spec=['t', 'e', 'ic'],
            options={})
    # selecting an isotherm at 10 eV
    Xarr = tab.Pt_DT['R_Array'][2:-1]*1.1 
    Yarr = np.array([10,  100, 1000])*11640
    X, Y = np.meshgrid(Xarr, Yarr, indexing='ij')

def test_second_derivatives():
    # check that Fxy == Fxy (assuming function C²)
    F0 = tab.Pt_DT.dFxy(X,Y)
    F1 = tab.Pt_DT.dFyx(X,Y)
    np.testing.assert_allclose(F0, F1, rtol=5e-2)

def test_derived_quantities_zero_order():
    """Checking all derived quantities can be computed (e.g. that dependencies are ok)"""
    for qtty_name, qtty_dict in eos.derived_quantities.iteritems():
        ctab = eos.EosMaterial(material, [qtty_name], spec=qtty_dict['spec'])
        for spec in qtty_dict['spec']:
            yield ctab.q[qtty_name, spec], X, Y

def test_sound_speed():
    """Checking that the isentropic sound speed computed by 2 methods is the same"""
    F0 = tab.q['Cs2'](X,Y)
    F1 = tab.q['Cs2_check'](X,Y)
    np.testing.assert_allclose(F0, F1)

def test_Cs2_larger_Ct2():
    """Check that isentropic sound speed square is larger then the isotermal one"""
    F0 = tab.q['Ct2'](X,Y)
    F1 = tab.q['Cs2'](X,Y)
    np.testing.assert_array_less(F0,F1)

def test_derived_quantities():
    """Check the thermodynamic relation:  Cv/Cp = Ks/Kt"""
    XYf = X,Y
    F0 = tab.q['Ks'](*XYf)/tab.q['Kt'](*XYf)
    F1 = tab.q['Cv'](*XYf)/tab.q['Cp'](*XYf)
    yield np.testing.assert_allclose, F0, F1

@nose.tools.nottest
def test_derived_quantities2():
    """Check the thermodynamic relation: 
                                      2
                C_v         T⋅V⋅α_exp
                ───  =  1 - ─────────
                C_p          C_p⋅K_t
        """
    XYf = X,Y
    F0 = tab.q['Ks'](*XYf)/tab.q['Kt'](*XYf)
    F1 = 1. - tab.q['alpha_exp'](*XYf)**2*Y/(X*tab.q['Cp'](*XYf)*tab.q['Kt'](*XYf))
    yield np.testing.assert_allclose, F0, F1

def test_thermal_expansion_coeff():
    """Checking the thermal expansion coefficient by two methods"""
    F0 = tab.q['alpha_exp_check'](X,Y)
    F1 = tab.q['alpha_exp'](X,Y)
    yield np.testing.assert_allclose, F0, F1, 1e-5, 1e-8


@nose.tools.nottest
def test_gruneisen_coeff():
    """Checking the Gruneisen coefficient by two methods"""
    yield np.testing.assert_allclose, tab.q['gruneisen_coef'](X,Y), tab.q['Gamma'](X,Y)

@nose.tools.nottest
def test_gammas():
    """Checking: Cp/Cv = γ·g/(γ·g-Γ²)"""
    d = {}
    for key in ['g', 'Gamma', 'gamma_3']:
        d[key] = tab.q[key](X,Y)
    F0 = tab.q['Cp'](X,Y)/tab.q['Cv'](X,Y)
    F1 = d['gamma_3']*d['g']/(d['gamma_3']*d['g'] - d['Gamma'])
    yield np.testing.assert_allclose, F0, F1, 1e-5



