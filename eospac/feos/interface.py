#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re

from .libpyfeos import  _get_calc_limits, _initialize, _init_mat,\
       _get_eos_all,_get_eos, _get_mat_par, _get_crit_point,\
       _get_soft_shere_par, _get_energy_offsets, _delete_material


from ..base import _pull_tables, TableBase, MaterialBase, GridBase
from copy import deepcopy


avalable_tabs = np.unique([tab.format(s=spec) \
                for tab in ['P{s}_DT', 'U{s}_DT', 'S{s}_DT', 'A{s}_DT', 'Zfc_DT',
                            'P{s}_DU{s}', 'T_DU{s}', 'S{s}_DU{s}', 'A{s}_DU{s}', 'Zfc_DU{s}']\
                for spec in ['t', 'iz', 'ec', 'e', 'ic']])

# define a global global_table_handles variable
if 'global_table_handles' not in globals():
    global_table_handles = None

class FeosTable(TableBase):
    _original_units = 'feos'

    def __init__(self, _name, table_handle=None, options={},info={}, units='cgs'):
        self._id = table_handle
        self._name = _name
        self._Nelements  = info['Nelements']
        self.update(info)
        super(FeosTable, self).__init__(_name, table_handle, options, units)
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
                res = _get_eos_all(self._id, self._Nelements,
                                          self.options['use_maxwell'], Xvals, Yvals)
                fVals_f = res[Fvar]*self._F_convert

                fVals = fVals_f.reshape(init_shape)
                return fVals
            else:
                raise NotImplemented
        else:
            raise NotImplemented

    def __getitem__(self, key):
        """ Overwriting default dict's __getitem__ method """
        if key in ['abar', "Mean_Atomic_Mass"]:
            return self['Atot']/self['Xtot']
        if key in ['zbar', "Mean_Atomic_Num"]:
            return self['Ztot']/self['Xtot']
        if key == "Normal_Density": key = "rho_ref"
        if key == "Modulus": key = "bulk_mod_ref"
        if key == "R_Array": key = "D_Array"
        if key == "Exchange_Coeff":
            return 0.0
        return super(FeosTable, self).__getitem__(key)

class FeosMaterial(MaterialBase):

    _default_options = {'use_maxwell': True, 'use_softspheres': True,
            'maxwell_temp_arr': None, 'max_materials': 1, 'debug': False,
            'rho_grid': None, 'temp_grid': None,
            'grid_subsample_default': 0, 'grid_kind': 'solid',
            'interpolation_mode_inversed': 'linear', 'precalculate': True}
    _backend = 'feos'
    _original_units = 'feos'

    def __init__(self, material=None, tables=['Pt_DT', 'Ut_DT'],
        options={'use_maxwell': True, 'use_softspheres': True,
            'maxwell_temp_arr': None, 'max_materials': 1, 'debug': False},
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
        self.options = self._validate_options(options)
        opt = self.options

        self.tables = _pull_tables(tables, spec, avalable_tabs)
        self.material = int(material)

        # making a global array to store table handles
        # This makes sure that FEOS is only initalized once
        global global_table_handles
        if global_table_handles is None:
            self._id = 1
            _initialize(opt['max_materials'], int(opt['debug']))
            global_table_handles = np.zeros(opt['max_materials'])
        else:
            for idx0, mat in enumerate(global_table_handles):
                if not mat:
                    self._id = idx0 + 1
                    break
            else:
                raise ValueError('Could not allocate {0} material.\n\
                    Please reinitialize all materials and increase the "max_materials" option!'.format(
                    self.material))

        global_table_handles[self._id-1] = self.material



        self.info = self._get_info_init()

        # default grid
        default_grid = GridBase(self.info['rho_ref'], kind=opt['grid_kind'],
                subsample_temp=opt['grid_subsample_default'],
                subsample_rho=opt['grid_subsample_default'])

        grid_units = self._set_units(units, 'Pt_DT')
        if opt['rho_grid'] is None:
            rho_grid = default_grid.rho_grid*grid_units.o2r('D', 'cgs', units)
        else:
            rho_grid = opt['rho_grid']

        if opt['temp_grid'] is None:
            temp_grid = default_grid.temp_grid*grid_units.o2r('T', 'cgs', units)
        else:
            temp_grid = opt['temp_grid']

        self.info['T_Array'] = temp_grid
        self.info['D_Array'] = rho_grid
        temp_grid_feos = temp_grid*grid_units.o2r('T', units, 'feos')
        rho_grid_feos = rho_grid*grid_units.o2r('D', units, 'feos')

        self.info.update(self._get_info_final(temp_grid_feos))
        if opt['precalculate']:
            self.precalculated = self._precalculate_on_grid(rho_grid_feos,temp_grid_feos)
            self.precalculated['D_Array'] = rho_grid_feos
            self.precalculated['T_Array'] = temp_grid_feos
        else:
            self.precalculated = {}
        for tab_idx, tab_key in enumerate(self.tables):
           this_info = deepcopy(self.info)
           if self.precalculated:
               F_Array = self.precalculated[tab_key.split('_')[0]]*\
                        grid_units.o2r(tab_key[0], 'feos', units)
               this_info['F_Array'] =  F_Array.T

           setattr(self, tab_key,
                FeosTable(tab_key,
                            self._id,
                            options=opt,
                            info=this_info,
                            units=units))
        self._init_base()

    def _precalculate_on_grid(self, rho_grid, temp_grid):
        R, T = np.meshgrid(rho_grid, temp_grid, indexing='ij')
        R_flat, T_flat = R.flatten(), T.flatten()
        #res = self._Nelements
        res_tmp = _get_eos_all(self._id, self._Nelements,
                self.options['use_maxwell'], R_flat, T_flat)
        res = {}

        for key in res_tmp:
            res[key] = res_tmp[key].reshape(R.shape)
        return res



    def _get_info_init(self):
        """ We have to preinitialize the material just to get the rho_ref """
        opt = self.options
        # preinitializint without Maxwell constructions to go faster
        Nelements, N_maxwell_iso, Maxwell_temp_reliable = _init_mat(self._id,
                self.material, False, False, 
                np.logspace(-4, 1, 100))

        info  = {}
        info.update(_get_mat_par(self._id, Nelements))
        _delete_material(self._id)

        return info

    def _get_info_final(self, temp_grid):
        """ We have to preinitialize the material just to get the density """
        opt = self.options
        T_mask = temp_grid<100. # eV 

        Nelements, N_maxwell_iso, Maxwell_temp_reliable = _init_mat(self._id,
                self.material, opt['use_maxwell'], opt['use_softspheres'],
                temp_grid)

        info  = {}
        self._Nelements = Nelements
        info['Nelements'] = Nelements
        info['N_maxwell_iso'] = N_maxwell_iso
        info['Maxwell_temp_reliable'] = Maxwell_temp_reliable
        info.update(_get_mat_par(self._id, Nelements))
        info.update(_get_energy_offsets(self._id))
        if opt['use_maxwell']:
            info.update(_get_crit_point(self._id))
            if opt['use_softspheres']:
                info.update(_get_soft_shere_par(self._id))

        return info

