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

    def save(self, filename, ext='sesame_bin', matid=None, energy_offset={}):
        """ Write the current material to a file.
        Parameters:
        -----------
          - filename: filename
          - ext: extension: sesame_bin, ionmix
          - energy_offset: dict with keys 'e', 'iz'
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
                   if spec in energy_offset:
                       U_offset = energy_offset[spec]
                   elif spec == 't' and 'e' in energy_offset and 'iz' in energy_offset:
                       U_offset = energy_offset['e'] + energy_offset['iz']
                   else:
                       U_offset = 0.
                   tabs[tab_id] = { 'rho': P_DT['R_Array']*units.o2r('D'),
                                    'temp': P_DT['T_Array']*units.o2r('T'),
                                    'U': U_DT['F_Array']*units.o2r('U') + U_offset,
                                    'P': P_DT['F_Array']*units.o2r('P'),
                                    'A': P_DT['F_Array']*units.o2r('A')}
                else:
                    print "Ignored {s} specie, as it doens't seem to be present in the table!".format(s=spec)
            # writing ionization
            if hasattr(self, 'Zfc_DT'):
                Zf = self.Zfc_DT
                tabs[601] = { 'rho': Zf['R_Array']*units.o2r('D'),
                              'temp': Zf['T_Array']*units.o2r('T'),
                              'Z':  Zf['F_Array']}
            else:
                print 'Ionization table not present!'
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



def _arithmetic_mean(vect):
    return 0.5*(vect[:-1] + vect[1:])

def _geometric_mean(vect):
    return (vect[:-1]*vect[1:])**0.5



class GridBase(object):
    # Normalized grid  adapted from table 3720
    # temperature grid  in eV
    rho_grid_init = np.concatenate((np.array([0]),
      np.logspace(-10,-7,15)[:-1], # 5 pt per decade
      np.logspace(-6,-1,50),       # 10 pt per decade
      np.array(
      [ 1.197e-01,  1.437e-01,  1.676e-01,  1.916e-01,  2.394e-01,  2.873e-01,
        3.352e-01,  3.830e-01,  4.307e-01,  4.789e-01,  5.267e-01,  5.744e-01,
        6.226e-01,  6.704e-01,  7.181e-01,  7.659e-01,  7.900e-01,  8.141e-01,
        8.378e-01,  8.619e-01,  8.859e-01,  9.096e-01,  9.289e-01,  9.481e-01,
        9.578e-01,  9.674e-01,  9.863e-01,  1.006e+00,  1.030e+00,  1.053e+00,
        1.101e+00,  1.149e+00,  1.245e+00,  1.341e+00,  1.437e+00,  1.532e+00,
        1.628e+00,  1.724e+00,  1.820e+00,  1.916e+00,  2.155e+00,  2.394e+00,
        2.634e+00,  2.873e+00,  3.113e+00,  3.352e+00,  3.592e+00,  3.830e+00,
        4.307e+00,  4.789e+00,  5.267e+00,  5.744e+00,  6.226e+00,  6.704e+00,
        7.181e+00,  7.659e+00,  8.141e+00,  8.619e+00,  9.096e+00,  9.578e+00,
        1.053e+01,  1.149e+01,  1.245e+01,  1.437e+01,  1.676e+01,  1.916e+01,
        2.155e+01,  2.394e+01,  2.873e+01,  3.352e+01,  3.830e+01,  4.789e+01,
        5.744e+01,  6.704e+01,  7.659e+01,  8.619e+01,  9.578e+01,  1.197e+02,
        1.437e+02,  1.916e+02,  2.873e+02,  3.830e+02,  5.744e+02,  7.659e+02,
        9.578e+02,  1.916e+03,  4.789e+03])))
    temp_grid_init = np.concatenate((np.array(
      [ 0.        ,  0.0062311 ,  0.01245704,  0.01868557,  0.02560997,
        0.0299055 ,  0.03488832,  0.03987113,  0.04486254,  0.04984536,
        0.05482818,  0.059811  ,  0.06480241,  0.06978522,  0.07476804,
        0.07975086,  0.08474227,  0.08969072,  0.09965636,  0.10962199,
        0.11958763,  0.13453608,  0.14948454,  0.16443299,  0.17938144,
        0.19931271,  0.22431271,  0.2492268 ,  0.29905498,  0.34888316,
        0.39871134,  0.44862543,  0.49845361,  0.54828179,  0.59810997,
        0.69785223,  0.79750859,  0.89690722,  0.99656357,  1.09621993,
        1.19587629,  1.34536082,  1.49484536,  1.6443299 ,  1.79381443,
        1.99312715,  2.24312715,  2.49226804,  2.74140893]),
        np.logspace(np.log10(3),1,15), np.arange(11,50,1), np.arange(50,80,2.5),
        np.arange(80,100,5), np.logspace(2,4,50),np.logspace(4,6,20)[1:] ))
    #temp_grid_init = np.logspace(-3, 6, 46)


    def __init__(self, rho_ref=None, kind='solid', subsample_rho=0, subsample_temp=0, temp_factor=1.0):
        self.temp_grid = self._subsample(
                                self.temp_grid_init*eV2K_cst*temp_factor, subsample_temp)
        if kind=='solid':
            if rho_ref is None:
                raise ValueError
            self.rho_grid = self._subsample(
                                self.rho_grid_init*rho_ref, subsample_rho)
        elif kind=='gas':
            rho_grid_init = np.concatenate( (np.array([0]),
                   np.logspace(-10,-7,15)[:-1],
                   np.logspace(-7, 1, 80)[:-1],
                   np.logspace(1, 2, 7)))
            self.rho_grid =  self._subsample(rho_grid_init, subsample_rho)
        elif kind=='log':
            self.temp_grid = np.logspace(-3, 6, 89)*eV2K_cst
            self.rho_grid = np.logspace(-8, 3.5, 46)

        else:
            raise NotImplemented

    @classmethod
    def _subsample(cls, grid, iterations=0, mean='geometric'):
        if iterations == 0:
            return grid
        if mean == 'arithmetic':
            grid_extra = _arithmetic_mean(grid)
        elif mean == 'geometric':
            grid_extra = _geometric_mean(grid)
        else:
            raise ValueError

        grid_extra = 0.5*(grid[:-1] + grid[1:])
        grid_new = np.concatenate((grid, grid_extra))
        grid_new = np.sort(grid_new)
        if iterations<=1:
            return grid_new
        else:
            return cls._subsample_arithmetic_mean(grid_new, iterations-1)



