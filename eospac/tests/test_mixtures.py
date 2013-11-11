#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import eospac.eospac as eos


def setup():
    global materials
    global tables_list
    global mix
    materials = [3720, 7593]
    tables_list = ['Pt_DT', 'Ut_DT']
    mix = eos.EospacMixture(materials, tables=tables_list)

def test_mixture1el():
    """Check that 1 element mixture gives the same result as single material"""
    mix1el = eos.EospacMixture(materials, ['Pt_DT'])
    mat_i =  mix1el.m3720
    mat_d = eos.EospacMaterial(3720, ['Pt_DT'])
    yield np.testing.assert_allclose, mat_i.Pt_DT['R_Array'], mat_d.Pt_DT['R_Array']
    yield np.testing.assert_allclose, mat_i.Pt_DT['F_Array'], mat_d.Pt_DT['F_Array']

def test_mixture2el_1dominant():
    """Check that a mixture of 2 materials with one in zero quantity is the same as single material"""

    Rmax =  mix.Pt_DT['Rmax'].min()
    Rmin =  mix.Pt_DT['Rmin'].max()
    Tmax =  mix.Pt_DT['Tmax'].min()
    Tmin =  mix.Pt_DT['Tmin'].max()
    #Xarr = np.logspace(np.log10(Rmin)+1, np.log10(Rmax)-1, 50)
    #Yarr = np.logspace(np.log10(Tmin)+1, np.log10(Tmax)-1, 80)
    Xarr = np.array([0.1])
    Yarr = np.array([7*11640])
    X, Y = np.meshgrid(Xarr, Yarr, indexing='ij')
    frac_dens1 = np.ones(X.shape)
    frac = np.array([0*frac_dens1,1.0*frac_dens1])
    #print mix.Pt_DT(X,Y,frac)[0]




