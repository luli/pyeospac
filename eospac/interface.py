#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re

from .libpyeospac import _create_tables, _load_tables, _interpolate,\
            _get_table_info, _set_option, _mix

import constants as cst

from .quantities import derived_quantities, DerivedQuantity

SMALL_H=1.0e-3

def _pull_tables(keys, spec):
    """Pull all necessary tables
    Parameters
    ----------
      - keys: list: a list containing table names including regular
        expressions, or derived quantities names. The '{s}' pattern will be
        expanded according to provided species
      - spec: list of species. Should be in ['t','e', 'ic']

    Returns
    -------
     - tabs: a list of tables in eospac.tables
    """

    keys_found = np.array(np.zeros(len(keys)), dtype='bool')
    keys = np.array(keys, dtype='object')

    requested_tables = []
    for key_idx, key in enumerate(keys):
       for quantity_name in derived_quantities:
           if re.match(key, quantity_name):
               requested_tables += derived_quantities.get_dependencies(quantity_name, spec)
               keys_found[key_idx] = True
       # if a derived quantity was matched for this key, go to the next
       # one
       if keys_found[key_idx]:
           continue

       for tab_name in cst.tables:
           for spec_el in spec:
               if re.match(key.format(s=spec_el), tab_name):
                   requested_tables.append(tab_name)
                   keys_found[key_idx] = True
    if not np.all(keys_found):
        raise ValueError("Keys {0} didn't match any table or derived quantity!".format(keys[~keys_found]))

    return list(set(requested_tables))  # making the list unique


class EosTable(dict):
    def __init__(self, table_handle, _name, info={}, options={}):
        self._id = table_handle
        self._name = _name
        for key in info:
            self[key] = info[key]
        self.options = options
        self._apply_options()

        # Try to set label:
        label = []
        label_dict = {'P': 'pressure', 'U':'internal energy', 'S': 'entropy',
                'A': 'free energy'}
        spec_dict = {'t': 'Total', 'e': 'Electron', 'i': 'Ion', 'ic': 'Ion+cc'}
        unit_dict = {'P': 'GPa', 'S': 'MJ.kg$^{-1}$.K$^{-1}$', 'U': 'MJ.kg$^{-1}$', 
                'A': 'MJ.kg$^{-1}$'}
        variable = _name[0]
        spec = _name.split('_')[0][1:]
        if spec in spec_dict:
            label.append(spec_dict[spec])
        else:
            label.append('??')
        if variable in label_dict:
            label.append(label_dict[variable])
            self.units = unit_dict[variable]
        else:
            label.append('??')
            self.units = '??'
        self.label = ' '.join(label)


    def _interpolate(self, X, Y, kind):
        if X.shape != Y.shape:
            raise ValueError('X and Y arguments should be ndarrays of the same shape!')

        init_shape = X.shape

        Xvals = np.array(X.flatten(), dtype='float64')
        Yvals = np.array(Y.flatten(), dtype='float64')
        fVals_f, dFx_f, dFy_f = _interpolate(self._id, Xvals, Yvals)

        fVals = fVals_f.reshape(init_shape)
        dFx = dFx_f.reshape(init_shape)
        dFy = dFy_f.reshape(init_shape)
        if kind == 'F':
            return fVals
        elif kind == 'dFx':
            return dFx
        elif kind == 'dFy':
            return dFy

        # computing the second derivative with finite differences
        Xl, Xr = X*(1-SMALL_H), X*(1+SMALL_H)
        Xl = np.fmax(Xl, 0)
        dX = Xr - Xl

        Yl, Yr = Y*(1-SMALL_H), Y*(1+SMALL_H)
        Yl = np.fmax(Yl, 0)
        dY = Yr - Yl

        if kind == 'dFxx':
            return (self._interpolate(Xr, Y, 'dFx') -\
                    self._interpolate(Xl, Y, 'dFx'))/dX
        elif kind == 'dFyy':
            return (self._interpolate(X, Yr, 'dFy') -\
                        self._interpolate(X, Yl, 'dFy'))/dY
        elif kind == 'dFxy':
            return (self._interpolate(Xr, Y, 'dFy') - \
                    self._interpolate(Xl, Y, 'dFy'))/dX
        elif kind == 'dFyx':
            return (self._interpolate(X, Yr, 'dFx') - \
                    self._interpolate(X, Yl, 'dFx'))/dY
        else:
            raise NotImplementedError

    def __call__(self, X, Y):
        return self._interpolate(X, Y, 'F')

    def dFx(self, X, Y):
        return self._interpolate(X, Y, 'dFx')

    def dFy(self, X, Y):
        return self._interpolate(X, Y, 'dFy')

    def dFxx(self, X, Y):
        return self._interpolate(X, Y, 'dFxx')

    def dFyy(self, X, Y):
        return self._interpolate(X, Y, 'dFyy')

    def dFxy(self, X, Y):
        return self._interpolate(X, Y, 'dFxy')

    def dFyx(self, X, Y):
        return self._interpolate(X, Y, 'dFyx')

    def __getitem__(self, key):
        """ Overwriting default dict's __getitem__ method """
        if key == "Z": key = "Mean_Atomic_Num"
        if key == "A": key = "Mean_Atomic_Mass"

        if key == "R_Array":
            info_keys = [cst.info[key]]*self['NR']
        elif key == "T_Array":
            info_keys = [cst.info[key]]*self['NT']
        elif key == "F_Array":
            info_keys = [cst.info[key]]*self['NR']*self['NT']
        else:
            info_keys = [cst.info[key]]

        info_val = _get_table_info(self._id, np.array(info_keys, dtype='int32'))
        if key == "F_Array":
            return info_val.reshape((self['NR'], self['NT']), order='F')
        elif len(info_val)>1:
            return info_val
        else:
            info_val = info_val[0]
            if key in ['NR', 'NT', 'Material_ID', 'Table_Type']:
                return int(info_val)
            else:
                return info_val

    def __repr__(self):
        np.set_printoptions(precision=3, threshold=40, edgeitems=2)
        out = ["="*80, 
               " "*30+"Table summary: {0}".format(self._name),
               "="*80,
               '\nParameters',
               '-'*15]
        for key in self.keys():
            try:
                out.append('{0:20} : {1}'.format(key, self[key]))
            except:
                pass
        out += ['\n'
                'Options',
                '-'*15]
        for key, val in self.options.iteritems():
            out.append('{0:20} : {1}'.format(key, val))

        out += ['='*80]
        np.set_printoptions() # going back to default options
        return '\n'.join(out)

    def max(self):
        return self['Fmax']

    def min(self):
        return self['Fmin']

    def _apply_options(self):
        """
        By default apply options in the same order as defined in
        eos_Interface.h
        """
        options_keys = sorted(self.options.keys(),
                                key=lambda k: cst.options[k])
        for key in options_keys:
            if self.options[key] is True: option_val = 1.0
            if self.options[key]:
                option_val = self.options[key]
            else:
                option_val = 0.0

            _set_option(self._id, cst.options[key],  option_val)
        _load_tables(np.array([self._id], dtype='int32'))



class EosMaterial:
    def __init__(self, material, tables=['Pt_DT', 'Ut_DT'], options={},
            spec=['t'], table_handles=None):
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
       self.material = int(material)

       self.tables = _pull_tables(tables, spec)
       self.options = options
       if table_handles is None:
           self._id_list = _create_tables(
                   np.array([cst.tables[key] for key in self.tables], dtype='int32'),
                   np.array([self.material]*len(self.tables), dtype='int32'))

       else:
           self._id_list = table_handles

       for tab_idx, tab_key in enumerate(self.tables):
           options = {}
           for opt_key, opt_vals in self.options.iteritems():
               if re.match(opt_key, tab_key):
                   options.update(opt_vals)
           setattr(self, tab_key,
                EosTable(self._id_list[tab_idx],
                            tab_key,
                            info=cst.info,
                            options=options))
       self.derived_quantity = DerivedQuantity(self)
       self.q = self.derived_quantity

    def get_table(self, table_name, spec=None):
        if spec is not None:
            return getattr(self, table_name.format(s=spec))
        else:
            return getattr(self, table_name)



class _EosMixtureSubclass:
    def __init__(self, table_handles):
        self._id_list = table_handles.copy()

    def _interpolate(self, X, Y, frac, kind):
        init_shape = X.shape
        init_frac_shape = (len(self._id_list), ) + init_shape
        if X.shape != Y.shape:
            raise ValueError('X and Y arguments should be ndarrays of the same shape!')
        if  init_frac_shape  !=  frac.shape:
            raise ValueError("Wrong shape for 'frac' argument: expected {0}, recieved {1}".format(
                init_frac_shape, frac.shape))

        Xvals = np.array(X.flatten(), dtype='float64')
        Yvals = np.array(Y.flatten(), dtype='float64')
        fracVals = np.array(frac.reshape((len(self._id_list), -1)), dtype='float64')
        fVals_f, dFx_f, dFy_f = _mix(self._id_list, Xvals, Yvals, fracVals)

        fVals = fVals_f.reshape(init_shape)
        dFx = dFx_f.reshape(init_shape)
        dFy = dFy_f.reshape(init_shape)
        if kind == 'F':
            return fVals
        elif kind == 'dFx':
            return dFx
        elif kind == 'dFy':
            return dFy

    def __call__(self, X, Y, frac):
        return self._interpolate(X, Y, frac, 'F')

    def dFx(self, X, Y, frac):
        return self._interpolate(X, Y, frac, 'dFx')
    
    def dFy(self, X, Y, frac):
        return self._interpolate(X, Y, frac, 'dFy')

    def __getitem__(self, key):
        if key == "Z": key = "Mean_Atomic_Num"
        if key == "A": key = "Mean_Atomic_Mass"

        info_key = cst.info[key]
        if key not in ['R_Array', 'T_Array', 'F_Array']:
            out = [_get_table_info(_id, np.array([info_key], dtype='int32'))[0]\
                    for _id in self._id_list]
            if key in ['NR', 'NT', 'Material_ID', 'Table_Type']:
                dtype='int32'
            else:
                dtype='float64'
            return np.array(out, dtype=dtype)
        else:
            raise NotImplementedError




class EosMixture:
    def __init__(self, materials, tables=['Pt_DT', 'Ut_DT'], options={}, spec=['t']):
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
        self.materials = np.array([int(mat) for mat in materials], dtype='int32')

        self.tables = _pull_tables(tables, spec)

        self.options = options
        Mat, Tabs = np.meshgrid(self.materials,
                     np.array([cst.tables[key]  for key in  self.tables]),
                     indexing='ij')
        Mat, Tabs = np.array(Mat, dtype='int32'), np.array(Tabs, dtype='int32')
        _id_list = _create_tables(Tabs.flatten(), Mat.flatten())
        self._id_list = _id_list.reshape(Mat.shape)

        for material_idx, material in enumerate(self.materials):
            material_options = {}
            for regexp_tuple, option_vals in options.iteritems():
                if type(regexp_tuple) is tuple:
                    mat_regexp, table_regexp = regexp_tuple
                    if re.match(str(mat_regexp), str(material)):
                        material_options[table_regexp] = option_vals
                else:
                    table_regexp = regexp_tuple
                    material_options[table_regexp] = option_vals

            setattr(self, 'm{0}'.format(material),
                   EosMaterial(material, tables=tables,
                       options= material_options,
                       table_handles=self._id_list[material_idx]
                       ))
        for table_idx, table_name in enumerate(self.tables):
            setattr(self, table_name, 
                    _EosMixtureSubclass(self._id_list[:,table_idx]))

        self.derived_quantity = DerivedQuantity(self)
        self.q = self.derived_quantity

    def get_table(self, table_name, spec=None):
        if spec is not None:
            return getattr(self, table_name.format(s=spec))
        else:
            return getattr(self, table_name)



