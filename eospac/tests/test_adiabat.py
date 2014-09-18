#!/usr/bin/python
# -*- coding: utf-8 -*-
from eospac import eV2K_cst
from eospac.adiabat import Adiabat
import numpy as np
from numpy.testing import assert_allclose



def setup():
    global rh_options
    rh_options = dict(backend='eospac',
            eos_pars={'.*': {'linear': False}},
                  material=3720)

def test_compare_with_original():
    """Compare with original SESAME implementation by Kerley """
    f = np.loadtxt('hugoniot_data/compression_adiabat_3720.txt')
    state0 = {'rho': 2.7, 'temp': 300}
    res_ref = {}
    for idx, key in enumerate(['rho', 'pres', 'eint', 'temp']):
        res_ref[key] = f[:,idx]
    #res_ref['temp'] *= eV2K_cst
    res_ref['pres'] = res_ref['pres']*100 # to Mbar -> GPa

    res = Adiabat.solve(state0, res_ref['rho'], **rh_options)
    #print res
    for key in res_ref:
        #print key, np.max(np.abs(res_ref[key] - np.array(res[key]))/np.fmin(np.abs(res_ref[key]), np.abs(np.array(res[key]))))
        yield assert_allclose, res_ref[key], np.array(res[key]), 1.0






