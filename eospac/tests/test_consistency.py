#!/usr/bin/python
# -*- coding: utf-8 -*-
import eospac as eos
import numpy as np


def setup():
    global tables_list, material
    global sestab


    tables_list = ['Pt_DT', 'Ut_DT']
    material = 7593
    sestab = eos.EosMaterial(material,
            tables=tables_list,
            options={},
            units='eospac')

def test_interpolation():
    """Checking that the interpolation gives the values in F_array"""
    for table in tables_list:
        ctable = getattr(sestab, table)
        R_array = ctable['R_Array']
        T_array = ctable['T_Array']
        F_array = ctable['F_Array']

        R,T = np.meshgrid(R_array, T_array, indexing='ij')

        F_array_itp = ctable(R, T)
        #F_array_itp = F_array_itp.T.reshape(F_array.shape)
        yield np.testing.assert_allclose, F_array, F_array_itp, 1e-15

#print al_eos.Ut_DT(X,Y)
#for tab in all_tables:
#    print tab, getattr(al_eos, tab).options
#print al_eos.Ut_DT


