#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright CNRS 2012
# Roman Yurchak (LULI)
# This software is governed by the CeCILL-B license under French law and
# abiding by the rules of distribution of free software.

from scipy.constants import physical_constants

from ..base import MaterialBase, TableBase, _pull_tables

R_CST = physical_constants['molar gas constant'][0]*1e7 # erg.K⁻¹.mol⁻¹

def Pt_DT(self, rho, temp):
    return rho*R_CST*temp/self['abar']

def Ut_DT(self, rho, temp):
    pres = rho*R_CST*temp/self['abar']
    eint = pres/((self['game']-1)*rho)
    return eint

def Pt_DUt(self, rho, eint):
    return (self['game'] - 1)*eint*rho

def T_DUt(self, rho, eint):
    pres = (self['game'] - 1)*eint*rho
    return pres*self['abar']/(rho*R_CST)

def Ut_DPt(self, rho, pres):
    eint = pres/((self['game']-1)*rho)
    return eint

def T_DPt(self, rho, pres):
    return pres*self['abar']/(rho*R_CST)

available_tables = {'Pt_DT': Pt_DT, 'Ut_DT': Ut_DT,
                    'Pt_DUt': Pt_DUt, 'T_DUt': T_DUt,
                    'Ut_DPt': Ut_DPt, 'T_DPt': T_DPt}


class GammaTable(TableBase):
    _original_units = 'cgs'
    def __init__(self, _name, table_handle=None,
            options={'game': None, 'gamc':None, 'abar':None}, units='cgs'):
        self.options = {}
        super(GammaTable, self).__init__(_name, table_handle, options, units)
        for key in ['game', 'gamc','abar']:
            self[key] = options[key] 

    def _interpolate(self, X, Y, kind):
        if kind == 'F':
            Xin, Yin = X*self._X_convert, Y*self._Y_convert
            return available_tables[self._name](self, Xin, Yin)*self._F_convert
        else:
            return self._differentiate(X, Y, kind)



class GammaMaterial(MaterialBase):
    def __init__(self, material=None, tables=None,
            options={'game': None, 'gamc': None, 'abar': None},
            spec=['t'], table_handles=None, units='cgs'):
        """
        Parameters:
        -----------
         - material: int: 4 digit SESAME material ID
         - tables: list: ['table1_id', 'table2_id', ..etc]
         - options: 
            game: defined in P = (γ_e - 1)ερ
            gamc: defined in γ_c = ρ/P * ∂P/∂ρ|_S
        """
        self.tables = _pull_tables(tables, spec, available_tables.keys())

        self.options = options
        if options['game']<1:
            raise ValueError('we should have game > 1')
        if options['gamc']<0:
            raise ValueError('we should have gamc > 0')
        self._init_base()

        for tab_idx, tab_key in enumerate(self.tables):
            setattr(self, tab_key,
                GammaTable(tab_key,
                            None,
                            options=options,
                            units=units))



