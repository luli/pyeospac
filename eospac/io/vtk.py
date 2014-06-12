#!/usr/bin/python
# -*- coding: utf-8 -*-
from interface import EosMaterial
from evtk.hl import gridToVTK
import numpy as np
from numpy.testing import assert_allclose

def _3d_reshape(arr):
    return arr.reshape(arr.shape+(1,))

def save_eostab_to_vtk(eosmat, filename):
    for view in ['DT']:
        for spec in ['e', 'i', 't']:
            pattern  = '_'.join([spec,view])
            ftables = [el for el in  eosmat.tables if pattern in el]
            if len(ftables):
                xl = []
                yl = []
                data = {}
                for tab_name in ftables:
                    tab = getattr(eosmat, tab_name)
                    xl.append(tab['R_Array'])
                    yl.append(tab['T_Array'])
                    data[tab_name] = _3d_reshape(tab['F_Array'])
                for x in xl:
                    assert_allclose(x,xl[0])
                for y in xl:
                    assert_allclose(y,xl[0])
                X, Y = np.meshgrid(x,y, indices='ij')
                X = _3d_reshape(X)
                Y = _3d_reshape(Y)
                Z = np.zeros(X.shape)

                #print X
                gridToVTK("./{0}".format(filename), X, Y, Z, pointData=data)

#def convert_to_vtk():
#    eosmat =  EosMaterial(3720, ['Pt_DT'])
#    save_eostab_to_vtk(eosmat, '3720')
