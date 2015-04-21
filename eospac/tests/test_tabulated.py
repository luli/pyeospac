#!/usr/bin/python
# -*- coding: utf-8 -*-
import eospac as eos
import numpy as np
import os.path
from numpy.testing import assert_allclose
from nose.plugins.skip import Skip, SkipTest

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

def test_sesascii():
    """ Conversion from ASCII to Binary for SESAME format """
    raise SkipTest

    matbase = 3719
    matid_new = int(9e4 + matbase)
    material = 'Fe'
    tab_options={'type': 'sesascii',
            'prescision': 'single',
            'path': '/home/rth/luli/eos/sesame/xsesame_ascii'}
    tables = ['Ut_DT', 'Pt_DT', 'At_DT']
  
    # Reading the ascii file 
    mat_ascii = eos.EosMaterial(matbase,
        tables=tables,
        options=tab_options,
        units='eospac', backend='tabulated')
    # Writing to binary
    filename = '/tmp/{0}_sesame_{1}.sesb'.format(material, matid_new)
    mat_ascii.save(filename, matid=matid_new)

    # Reopening the binary file
    mat_bin0 = eos.EosMaterial(matbase,
        tables=tables,
        options={},
        units='eospac', backend='eospac')

    mat_bin1 = eos.EosMaterial(matid_new,
        tables=tables,
        options={},
        units='eospac', backend='eospac')

    for tab_name in tables:
        tab0 =  getattr(mat_bin0, tab_name)
        tab1 =  getattr(mat_bin1, tab_name)

        for key in ['D_Array', 'T_Array', 'F_Array']:
            if not np.allclose(tab0[key], tab1[key]):
                print tab_name, key, 'failed'

            yield assert_allclose, tab0[key], tab1[key], 1e-5

#    yield assert_allclose, eosmat.Pt_DT(rho, temp), rho*R_CST*temp/abar
#    yield assert_allclose, eosmat.Pt_DT.dFx(rho, temp), R_CST*temp/abar
#    yield assert_allclose,   eosmat.Pt_DT.dFxx(rho, temp), np.zeros(Npt), 1e-8, 20
#    yield assert_allclose,   eosmat.Pt_DT.dFxy(rho, temp), R_CST/abar*np.ones(Npt)
#    yield assert_allclose, eosmat.Pt_DT.dFy(rho, temp), rho*R_CST/abar

#print al_eos.Ut_DT(X,Y)
#for tab in all_tables:
#    print tab, getattr(al_eos, tab).options
#print al_eos.Ut_DT


