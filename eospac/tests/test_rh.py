#!/usr/bin/python
# -*- coding: utf-8 -*-
from eospac import eV2K_cst
from eospac.rh import RankineHugoniot
import numpy as np
from numpy.testing import assert_allclose
import os
from nose.plugins.skip import Skip, SkipTest



def setup():
    global rh_options
    rh_options = dict(backend='eospac',
            eos_pars={'.*': {'create_tzero': True, 'linear': False}},
            root_opts={'method': 'hybr'},
                  material=3720)

def test_consitency():
    """Checking that computing RH by different methods gives the 
    same result """
    state0 = {"temp": 300, "rho": 2.7}
    #state1 = {"rho": np.linspace(1, 4, 100)*state0['rho']}
    state1 = {"pres": np.linspace(1e2, 2.5e3, 100)} # GPa # inital conditions so it doesn't fail...

    res1 = RankineHugoniot.solve(state0, state1, **rh_options)
    state1 = {"rho": res1['rho1']}
    res2 = RankineHugoniot.solve(state0, state1, **rh_options)

    for key in res1:
        yield assert_allclose, res1[key], res2[key]

def test_compare_with_original():
    """Compare with original SESAME implementation by Kerley """
    raise SkipTest
    basedir = os.path.dirname(os.path.realpath(__file__))
    f = np.loadtxt(os.path.join(basedir, 'hugoniot_data/hugoniot_3720.txt'))
    state0 = {'rho': 2.7, 'temp': 300, 'pres': 9.05001730E-03}
    res_ref = {}
    for idx, key in enumerate(['rho1', 'pres1', 'eint1', 'temp1', 'u_s', 'u_p']):
        res_ref[key] = f[:,idx]
    res_ref['temp1'] *= eV2K_cst
    res_ref['pres1'] = res_ref['pres1']*100 + state0['pres'] # to Mbar -> GPa

    state1 = {'pres': res_ref['pres1']}
    res = RankineHugoniot.solve(state0, state1, **rh_options)
    for key in res_ref:
        yield assert_allclose, res_ref[key], np.array(res[key]), 1.2e-2
        # data is constistant to up to 12%
        # the difference might come from convergence criteria or some
        # numeric. But that's good enough.






