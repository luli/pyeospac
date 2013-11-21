#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import physical_constants
import re

from .quantities import derived_quantities, DerivedQuantity

eV2K_cst = physical_constants['electron volt-kelvin relationship'][0]

def _pull_tables(keys, spec, valid_tables=[]):
    """Pull all necessary tables
    Parameters
    ----------
      - keys: list: a list containing table names including regular
        expressions, or derived quantities names. The '{s}' pattern will be
        expanded according to provided species
      - spec: list of species. Should be in ['t','e', 'ic', 'iz']
      - valid_tables: list of valid tables # cst.tables

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

       for tab_name in valid_tables:
           for spec_el in spec:
               if re.match(key.format(s=spec_el), tab_name):
                   requested_tables.append(tab_name)
                   keys_found[key_idx] = True
    if not np.all(keys_found):
        raise ValueError("Keys {0} didn't match any table or derived quantity!".format(keys[~keys_found]))

    return list(set(requested_tables))  # making the list unique

class EosUnits(dict):
    def __init__(self, original, requested):
        """
        Make the translation between different EoS units.
        For the moment the supported units are 'cgs', 'eospac', 'feos'
        Parameters:
        -----------
         - original: original unit sytem for the specific backend
         - requested: requested unit system
        """
        
        # Getting all necessary constants

        # Translation dicts in cgs units
        self['cgs']    = {'P': 1., 'T': 1., 'D': 1., 'U': 1.,
                        'A': 1., 'S': 1., 'Z': 1., 'V': 1.}

        self['feos']   = {'P': 1., 'T': eV2K_cst, 'D': 1., 'U': 1.,
                'A': 1., 'S': 1./eV2K_cst, 'Z': 1., 'V': 1.}

        self['eospac'] = {'P': 1.e10, 'T': 1., 'D': 1., 'U': 1.e10,
                'A': 1.e10, 'S': 1.e10, 'Z': 1., 'V': 1.}
        if original not in self:
            raise KeyError('Unit system {0} not implemented!'.format(original))
        if requested not in self:
            raise KeyError('Unit system {0} not implemented!'.format(requested))

        self.original = original
        self.requested = requested

    def o2r(self, var, original=None, requested=None):
        """ Original to requested unit system """
        if original is None:
            original = self.original
        if requested is None:
            requested = self.requested
        try:
            return self[original][var]/self[requested][var]
        except KeyError:
            print 'eospac/base.py:EosUnits.o2r Warning: variable {0} not defined within \n\
              EosUnits. You are either trying to access one of less used EOSPAC variables,\n\
              or something went wrong! Falling back to a converion factor of 1...'.format(var)
            return 1.
        else:
            raise

    def table_units(self, table_name):
        """
        Given a table name return necessary unit conversions
        """
        tab_vars = re.sub('[^A-Z]', '', table_name)
        if len(tab_vars) == 3:
            Fut, Xut, Yut = tab_vars
            return self.o2r(Fut), 1./self.o2r(Xut), 1./self.o2r(Yut)
        elif len(tab_vars) == 2:
            Fut, Xut = tab_vars
            return self.o2r(Fut), 1./self.o2r(Xut), None
        else:
            raise ValueError("Couldn't parse the table name!")



def _set_units(self, units, _name):
        self._requested_units = units
        physical_units = EosUnits(self._original_units, self._requested_units)
        self._F_convert, self._X_convert, self._Y_convert = physical_units.table_units(_name)
        return physical_units


direct_tables = ['P{s}_DT', 'U{s}_DT', 'A{s}_DT']



class TableBase(dict):
    _original_units = 'cgs'

    def __init__(self, _name, table_handle=None, options={}, units='cgs'):
        self._id = table_handle
        self._name = _name

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
        # setting units
        self._set_units(units, self._name)

    _set_units = _set_units


    #def _interpolate(self, X, Y, kind):
    #    if X.shape != Y.shape:
    #        raise ValueError('X and Y arguments should be ndarrays of the same shape!')

    #    init_shape = X.shape

    #    Xvals = np.array(X.flatten(), dtype='float64')
    #    Yvals = np.array(Y.flatten(), dtype='float64')
    #    fVals_f, dFx_f, dFy_f = _interpolate(self._id, Xvals, Yvals)

    #    fVals = fVals_f.reshape(init_shape)
    #    dFx = dFx_f.reshape(init_shape)
    #    dFy = dFy_f.reshape(init_shape)
    #    if kind == 'F':
    #        return fVals
    #    elif kind == 'dFx':
    #        return dFx
    #    elif kind == 'dFy':
    #        return dFy

    #    # computing the second derivative with finite differences
    #    Xl, Xr = X*(1-SMALL_H), X*(1+SMALL_H)
    #    Xl = np.fmax(Xl, 0)
    #    dX = Xr - Xl

    #    Yl, Yr = Y*(1-SMALL_H), Y*(1+SMALL_H)
    #    Yl = np.fmax(Yl, 0)
    #    dY = Yr - Yl

    #    if kind == 'dFxx':
    #        return (self._interpolate(Xr, Y, 'dFx') -\
    #                self._interpolate(Xl, Y, 'dFx'))/dX
    #    elif kind == 'dFyy':
    #        return (self._interpolate(X, Yr, 'dFy') -\
    #                    self._interpolate(X, Yl, 'dFy'))/dY
    #    elif kind == 'dFxy':
    #        return (self._interpolate(Xr, Y, 'dFy') - \
    #                self._interpolate(Xl, Y, 'dFy'))/dX
    #    elif kind == 'dFyx':
    #        return (self._interpolate(X, Yr, 'dFx') - \
    #                self._interpolate(X, Yl, 'dFx'))/dY
    #    else:
    #        raise NotImplementedError

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

    def _differentiate(self, X, Y, kind, drel=1e-1):
        itp = self.__call__
        Xl, Xr = X*(1-drel), X*(1+drel)
        Xl = np.fmax(Xl, 0)
        dX = 0.5*(Xr - Xl)

        Yl, Yr = Y*(1-drel), Y*(1+drel)
        Yl = np.fmax(Yl, 0)
        dY = 0.5*(Yr - Yl)

        if kind == 'dFx':
            return (itp(Xr, Y) - itp(Xl, Y))/(2*dX)
        elif kind == 'dFy':
            return (itp(X, Yr) - itp(X, Yl))/(2*dY)
        elif kind == 'dFxx':
            return (itp(Xr, Y)- 2*itp(X,Y) + itp(Xl, Y))/dX**2
        elif kind == 'dFyy':
            return (itp(X, Yr)- 2*itp(X,Y) + itp(X, Yl))/dY**2
        elif kind in ['dFxy', 'dFyx']:
            return (itp(Xr, Yr) - itp(Xr, Yl) - itp(Xl,Yr) + itp(Xl, Yl))/(4*dX*dY)
        else:
            raise NotImplementedError


    def __repr__(self):
        np.set_printoptions(precision=3, threshold=40, edgeitems=2)
        out = ["="*80, 
               " "*30+"Table summary: {0}".format(self._name),
               "="*80,
               '\nParameters',
               '-'*15]
        
        for key in sorted(self.keys()):
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


class MaterialBase(dict):
    def _init_base(self):
        """ Finish with general initialization. See backend specific __init__
        for more details (e.g. EospacMaterial.__init___)
        """

        self.derived_quantity = DerivedQuantity(self)
        self.q = self.derived_quantity

    def get_table(self, table_name, spec=None):
        if spec is not None:
            return getattr(self, table_name.format(s=spec))
        else:

            return getattr(self, table_name)

    def _get_state_DT(self, X, Y, spec='t'):
        rho = X
        temp = Y
        pres = self.get_table('P{s}_DT', spec)(X,Y)
        eint = self.get_table('U{s}_DT', spec)(X,Y)
        return {'eint': eint, 'rho': rho, 'temp': temp, 'pres': pres}

    def _get_state_DU(self, X, Y, spec='t'):
        rho = X
        eint = Y
        pres = self.get_table('P{s}_DU{s}', spec)(X,Y)
        temp = self.get_table('T_DU{s}', spec)(X,Y)
        return {'eint': eint, 'rho': rho, 'temp': temp, 'pres': pres}

    def _get_state_DP(self, X, Y, spec='t'):
        rho = X
        pres = Y
        eint = self.get_table('U{s}_DP{s}', spec)(X,Y)
        temp = self.get_table('T_DP{s}', spec)(X,Y)
        return {'eint': eint, 'rho': rho, 'temp': temp, 'pres': pres}

    def get_state(self, spec='t', mode=None, **args):
        """
        Get full EoS state.
        Mode should be None, DT, DU or DP. If it is None, an attempt is made to determine
        it based on provided args.
        Args keys can be only 'rho', 'temp', 'eint', 'pres'
        """
        mode_dict = {'DT': ('rho', 'temp'),
                     'DU': ('rho', 'eint'),
                     'DP': ('rho', 'pres'),
                    }
        if mode is None:
            # determining the mode automatically based on the provided arguments
            for key, req_vars in mode_dict.iteritems():
                if np.all([var in args for var in req_vars]):
                    mode = key
        if mode is None:
            raise ValueError('Failed to determine EoS mode, please provide the mode or give\
                              consistent set of input arguments.')
        req_vars = mode_dict[mode]
        X = args[req_vars[0]]
        Y = args[req_vars[1]]
        # trying to determine
        try:
            f =  getattr(self, '_get_state_{0}'.format(mode))
        except:
            raise KeyError("Mode {0} not implemented. Current valid modes are DT, DU.")
        return f(X, Y, spec)

    def _validate_options(self, options_in):
        """
        Check that provided options are consistent.
        """
        # check that given options don't contain unknown key
        for key in options_in:
            if key not in self._default_options:
                raise KeyError("Unknown option key {0}. Accepted options are: {1}".format(
                                    key, str(self._default_options.keys())))
        # setting to default version if key was not provided
        for key, default_val in self._default_options.iteritems():
            if key not in options_in:
                options_in[key] = default_val
        if self._backend == 'feos':
            pass
        return options_in

    _set_units = _set_units

    def save(self, filename, ext='sesame_bin', matid=None):
        """ Write the current material to a file.
        Parameters:
        -----------
          - filename: filename
          - ext: extension: sesame_bin, ionmix
        """

        if ext == 'sesame_bin':
            from .eospac.libsesio import _write_sesbin
            P_DT = self.get_table('P{s}_DT', spec='t')
            props = {}
            if matid is None:
                matid = P_DT['Material_ID']
            for key in ['Mean_Atomic_Num', 'Mean_Atomic_Mass',
                    'Normal_Density', 'Modulus', 'Exchange_Coeff']:
                props[key] = P_DT[key]
            units = EosUnits(self._requested_units, 'eospac')
            tabs = {}
            for tab_id, spec in [(301, 't'), (303, 'ic'),
                                 (304, 'e'), (305, 'iz')]:
                if np.all([key.format(s=spec) in self.tables for key in direct_tables]):
                   P_DT = self.get_table('P{s}_DT', spec=spec)
                   U_DT = self.get_table('U{s}_DT', spec=spec)
                   A_DT = self.get_table('A{s}_DT', spec=spec)
                   tabs[tab_id] = { 'rho': P_DT['R_Array']*units.o2r('D'),
                                    'temp': P_DT['T_Array']*units.o2r('T'),
                                    'U': U_DT['F_Array']*units.o2r('U'),
                                    'P':  P_DT['F_Array']*units.o2r('P'),
                                    'A': A_DT['F_Array']*units.o2r('A')}
                else:
                    print "Ignored {s} specie, as it doens't seem to be present in the table!".format(s=spec)
            # writing ionization
            Zf = self.get_table('Zfc_DT')
            tabs[601] = { 'rho': Zf['R_Array']*units.o2r('D'),
                          'temp': Zf['T_Array']*units.o2r('T'),
                          'Z':  Zf['F_Array']}
            #for tab_id in [601]:
            #    ctab = tabs[tab_id]
            #    ctab['rho'] = np.log10(np.fmax(1e-10, ctab['rho']))
            #    ctab['temp'] = np.log10(np.fmax(1e-10, ctab['temp']))

            _write_sesbin(filename, matid, props, tabs)
        else:
            raise NotImplemented

try:
    from .io import save2vtk
    MaterialBase.save2vtk = save2vtk
except ImportError:
    pass





class GridBase(object):
    # Normalized grid from 3720 table
    rho_grid_init = np.array([0.000e+00, 1e-10, 5.0e-10, 1.0e-9, 5.0e-9,
        1.0e-8, 5.0e-8, 1.0e-7, 5.0e-7,
        1.000e-6, 2.586e-06, 5.172e-06, 1.293e-05, 2.586e-05, 5.172e-05,
        1.293e-04, 2.586e-04, 5.172e-04, 1.293e-03, 2.586e-03, 3.879e-03, 6.465e-03,
        1.034e-02, 1.551e-02, 2.586e-02, 3.879e-02, 6.465e-02, 1.034e-01, 1.551e-01,
        2.068e-01, 2.586e-01, 3.232e-01, 3.879e-01, 4.525e-01, 5.172e-01, 6.465e-01,
        7.758e-01, 9.051e-01, 1.034e+00, 1.163e+00, 1.293e+00, 1.422e+00, 1.551e+00,
        1.681e+00, 1.810e+00, 1.939e+00, 2.068e+00, 2.133e+00, 2.198e+00, 2.262e+00,
        2.327e+00, 2.392e+00, 2.456e+00, 2.508e+00, 2.560e+00, 2.586e+00, 2.612e+00,
        2.663e+00, 2.715e+00, 2.780e+00, 2.844e+00, 2.974e+00, 3.103e+00, 3.362e+00,
        3.620e+00, 3.879e+00, 4.137e+00, 4.396e+00, 4.655e+00, 4.913e+00, 5.172e+00,
        5.818e+00, 6.465e+00, 7.112e+00, 7.758e+00, 8.405e+00, 9.051e+00, 9.698e+00,
        1.034e+01, 1.163e+01, 1.293e+01, 1.422e+01, 1.551e+01, 1.681e+01, 1.810e+01,
        1.939e+01, 2.068e+01, 2.198e+01, 2.327e+01, 2.456e+01, 2.586e+01, 2.844e+01,
        3.103e+01, 3.362e+01, 3.879e+01, 4.525e+01, 5.172e+01, 5.818e+01, 6.465e+01,
        7.758e+01, 9.051e+01, 1.034e+02, 1.293e+02, 1.551e+02, 1.810e+02, 2.068e+02,
        2.327e+02, 2.586e+02, 3.232e+02, 3.879e+02, 5.172e+02, 7.758e+02, 1.034e+03,
        1.551e+03, 2.068e+03, 2.586e+03, 5.172e+03, 1.293e+04, 2.586e+04, 5.172e+04])
    temp_grid_init = np.array([0.000e+00, 7.253e+01, 1.450e+02, 2.175e+02, 2.981e+02, 3.481e+02, 4.061e+02, 
        4.641e+02, 5.222e+02, 5.802e+02, 6.382e+02, 6.962e+02, 7.543e+02, 8.123e+02,
        8.703e+02, 9.283e+02, 9.864e+02, 1.044e+03, 1.160e+03, 1.276e+03, 1.392e+03,
        1.566e+03, 1.740e+03, 1.914e+03, 2.088e+03, 2.320e+03, 2.611e+03, 2.901e+03, 
        3.481e+03, 4.061e+03, 4.641e+03, 5.222e+03, 5.802e+03, 6.382e+03, 6.962e+03,
        8.123e+03, 9.283e+03, 1.044e+04, 1.160e+04, 1.276e+04, 1.392e+04, 1.566e+04,
        1.740e+04, 1.914e+04, 2.088e+04, 2.320e+04, 2.611e+04, 2.901e+04, 3.191e+04,
        3.481e+04, 4.061e+04, 4.641e+04, 5.802e+04, 7.543e+04, 9.283e+04, 1.160e+05,
        1.450e+05, 1.740e+05, 2.320e+05, 2.901e+05, 3.481e+05, 4.641e+05, 5.802e+05,
        9.283e+05, 1.160e+06, 1.740e+06, 2.901e+06, 4.641e+06, 6.962e+06, 9.283e+06,
        1.160e+07, 1.740e+07, 2.901e+07, 6.962e+07, 1.160e+08, 2.320e+08, 5.802e+08,
        1.160e+09, 2.320e+09, 5.802e+09, 1.160e+10])

    def __init__(self, rho_ref, kind='solid', subsample_rho=1, subsample_temp=1, temp_factor=1.0):
        if kind=='solid':
            self.rho_grid = self._subsample_arithmetic_mean(
                                self.rho_grid_init*rho_ref/2.7, subsample_rho)
            self.temp_grid = self._subsample_arithmetic_mean(
                                self.temp_grid_init*temp_factor, subsample_temp)
        else:
            raise NotImplemented

    @classmethod
    def _subsample_arithmetic_mean(cls, grid, iterations=1):
        if iterations == 0:
            return grid
        grid_extra = 0.5*(grid[:-1] + grid[1:])
        grid_new = np.concatenate((grid, grid_extra))
        grid_new = np.sort(grid_new)
        if iterations<=1:
            return grid_new
        else:
            return cls._subsample_arithmetic_mean(grid_new, iterations-1)



