#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re
from copy import deepcopy

from ..base import _pull_tables, TableBase, MaterialBase, GridBase


avalable_tabs = np.unique([tab.format(s=spec) \
                for tab in ['P{s}_DT', 'U{s}_DT', 'A{s}_DT', 'Zfc_DT']\
                for spec in ['t', 'ic', 'e', 'iz']])


supported_formats = {'ionmix': ('feos',)}

class TabulatedTable(TableBase):
    #_original_units = 'feos'

    def __init__(self, _name, table_handle=None, options={},info={}, units='cgs'):
        self._id = table_handle
        self._name = _name
        self.update(info)
        super(TabulatedTable, self).__init__(_name, table_handle, options, units)
        self.options = options

    def _interpolate(self, X, Y, kind):
        if type(X) in [float, int]:
           X = np.array(X)
           Y = np.array(Y)
        if X.shape != Y.shape:
            raise ValueError('X and Y arguments should be ndarrays of the same shape!')

        init_shape = X.shape

        Xvals = np.array(X.flatten(), dtype='float64')*self._X_convert
        Yvals = np.array(Y.flatten(), dtype='float64')*self._Y_convert
        Fvar, XY_vars =  self._name.split('_')
        if XY_vars == 'DT':
            if kind == 'F':
                return
            else:
                raise NotImplemented
        else:
            raise NotImplemented

    def __getitem__(self, key):
        """ Overwriting default dict's __getitem__ method """
        if key == "Mean_Atomic_Num" : key = "zbar"
        if key == "Mean_Atomic_Mass": key = "abar"
        if key == "Normal_Density": key = "rho_ref"
        if key == "R_Array": key = "D_Array"
        return super(TabulatedTable, self).__getitem__(key)

class TabulatedMaterial(MaterialBase):

    _default_options = {'path': None, 'type': 'ionmix', 'abar': -1., 'zbar': -1,
            'rho_ref': -1, 'Modulus': 0., 'Exchange_Coeff': 0.0, 'create_tzero': "zero"}
    _backend = 'tabulated'
    #_original_units = 'feo'

    def __init__(self, material=None, tables=['Pt_DT', 'Ut_DT'], options={},
            spec=['t'], units='cgs'):
        """
        Parameters:
        -----------
         - material: int: 4 digit SESAME material ID
         - tables: list: ['table1_id', 'table2_id', ..etc]
         - options: dict: {'tableid_regexpr': {'opt1': optval}, etc}
            For example {".t_DT": {'create_tzero': True}} would apply the
            EOS_CREATE_TZERO option to both 'Pt_DT' and 'Ut_DT' tables. The
            tableid_regexpr accepts regular expressions.
            [see  help(re) for more details].

        """
        self._id = -1
        options = self._validate_options(options)
        self.options = self._validate_options2(options)
        opt = self.options

        self.tables = [key for key in avalable_tabs]

        if material is not None:
            self.material = int(material)

        self.info = self._get_info()
        info = self.info
        ftype = opt['type']
        for key in ['type', 'path', 'zbar', 'abar', 'rho_ref',
                        'Modulus', 'Exchange_Coeff']:
            self.info[key] = opt[key]
            del opt[key]


        self._original_units = supported_formats[ftype][0]
        self._requested_units = units
        grid_units = self._set_units(units, 'Pt_DT')

        translation_dict = {'Pe_DT': 'pele', 'Pic_DT': 'pion', 'Ue_DT': 'eele',
                            'Uic_DT': 'eion', 'Zf_DT': 'zbar'}
        out = {}
        if ftype == 'ionmix':
            import opacplot2 as opp
            d = opp.OpacIonmix(info['path'], info['abar']/opp.NA, twot=True, man=True, verbose=False)
            D_array = d.dens*grid_units.o2r('D')
            T_array = d.temps*grid_units.o2r('T')

            for key in avalable_tabs:
                if 't' not in key:
                    if 'A' not in  key:
                        out[key] = getattr(d, translation_dict[key])*grid_units.o2r(key[0])
                    else:
                        out[key] = np.zeros(d.pele.shape)

            out['Pt_DT'] = out['Pe_DT'] + out['Pic_DT']
            out['Ut_DT'] = out['Ue_DT'] + out['Uic_DT']
            out['At_DT'] = out['Ae_DT'] + out['Aic_DT']

        for tab_idx, tab_key in enumerate(self.tables):
           this_info = deepcopy(self.info)
           this_info['D_Array'] = D_array
           self._create_tzero(D_array, T_array, out[tab_key])
           if opt['create_tzero']:
              this_info['T_Array'] = np.concatenate((np.zeros(1), T_array))
              this_info['F_Array'] = self._create_tzero(D_array, T_array, out[tab_key],
                                        method=opt['create_tzero'])
           else:
              this_info['T_Array'] = T_array
              this_info['F_Array'] = out[tab_key]

           setattr(self, tab_key,
                TabulatedTable(tab_key,
                            self._id,
                            options=opt,
                            info=this_info,
                            units=units))
        self._init_base()

    @staticmethod
    def _create_tzero(D_array, T_array, F_array, method='linear'):
        F_shape = F_array.shape
        F_array_new = np.ones((F_shape[0], F_shape[1]+1))
        F_array_new[:,1:] = F_array
        if method == 'linear':
            x = T_array[1:3] # first 2 non null points in temperature
            y = F_array_new[:, 1:3] 
            poly = np.polyfit(x, y.T, 1)
            F_array_new[:,0] = np.polyval(poly, 0) # evaluate @ 0
            return F_array_new
        elif method == 'zero':
            return F_array_new




    def _validate_options2(self, opt):
        if opt['type'] not in supported_formats:
            raise ValueError("Unvalid table type: {0}!".format(opt['type']))
        for key in ['zbar', 'abar', 'rho_ref']:
            if opt[key] < 0:
                raise ValueError('Parameter {0} should be specified!'.format(key))
        if opt['create_tzero'] is True:
            opt['create_tzero'] = 'zero'
        return opt

    def _get_info(self):
        info = {}
        info['Material_ID'] = self.material
        return info
