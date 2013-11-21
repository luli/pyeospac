#!/usr/bin/pythonm
# -*- coding: utf-8 -*-
import numpy as np

class DerivedQuantityContainter(dict):

    def add_quantity(self, **kwargs ):
        self[kwargs['name']] = {'name': kwargs['name'],
                      'dependencies': kwargs['dependencies'],
                      'function': kwargs['function'],
                      'args': kwargs['args'],
                      'description': kwargs['description'],
                      'spec': kwargs['spec']}
    def __repr__(self):
        out = ['List of derived quantities',
                '-'*30]
        for key, val in self.iteritems():
            out.append('   * {0:15} : {1:6} : {2}'.format(key, ' '.join(val['spec']),
                        val['description']))
        return '\n'.join(out)

    def get_dependencies(self, quantity, spec=['t']):
        """Get all the tables required to compute the given list of derived quantities."""
        #req_tabs = []
        #for key, val in self.iteritems():
        #    if key in quantities:
        #        req_tabs += val['dependencies']

        #req_tabs_all = []
        #for spec_el in spec:
        #    req_tabs_all += [el.format(s=spec_el) for el in req_tabs]

        #return list(set(req_tabs_all)) # making this list unique
        required_tables = []
        dependacies = self[quantity]['dependencies']
        for key in dependacies:
            if key in self:
                required_tables += self.get_dependencies(key, spec=spec)
            else:
                required_tables += [key.format(s=el) for el in spec]
        return list(set(required_tables))


derived_quantities = DerivedQuantityContainter()

class DerivedQuantity(dict):

    def __init__(self, tab):
        self.table = tab

    def __getitem__(self, args):

        if type(args) is tuple and len(args) == 2:
            key, spec = args
        elif type(args).__name__ in ['str']:
            key = args
            spec = 't'
        else:
            raise KeyError('Wrong number of keys, expected 1 or 2, recieved {0}'.format(len(args)))

        if key in derived_quantities:
            element = derived_quantities[key]
            #dependencies = [name.format(s=spec) for name in element['dependencies']]
            return lambda *XYt: element['function'](self.table, spec, *XYt)
        else:
            raise KeyError('Key {0} not in {1}'.format(key, str(self.keys())))


def add_quantity(**kwargs):
    def inner_decorator(function):
        if 'name' not in kwargs:
            kwargs['name'] = function.func_name
        if 'dependencies' not in kwargs:
            kwargs['dependencies'] = []
        if 'spec' not in kwargs:
            kwargs['spec'] = ['t', 'e', 'ic']
        kwargs['description'] = function.__doc__
        kwargs['function'] = function
        derived_quantities.add_quantity(**kwargs)
        return function
    return inner_decorator


# Isentropic quantities
# ---------------------
@add_quantity(name='Cs2_old', args='DT', dependencies=['P{s}_DS{s}', 'S{s}_DT'])
def isentropic_sound_speed2_old(tab, spec, *XYf):
    """Isentropic sound speed square"""
    XYf = list(XYf)
    XYf[1] = tab.get_table('S{s}_DT', spec)(*XYf)
    return tab.get_table('P{s}_DS{s}', spec).dFx(*XYf)

@add_quantity(name='Cs2', args='DT', dependencies=['P{s}_DT', 'P{s}_DU{s}', 'U{s}_DT'])
def isentropic_sound_speed2(tab, spec, *XYf):
    """Isentropic sound speed gamma: ρ★Cs²/P"""
    XYf_DT = list(XYf)
    dens = XYf[0]
    pres = tab.get_table('P{s}_DT', spec)(*XYf_DT)
    eint = tab.get_table('U{s}_DT', spec)(*XYf_DT)
    XYf_DU = XYf_DT[:]
    XYf_DU[1] = eint[:]

    dP_drho = tab.get_table('P{s}_DT', spec).dFx(*XYf_DT)
    dP_deint = tab.get_table('P{s}_DU{s}', spec).dFy(*XYf_DU)

    return  dP_drho+ pres*dP_deint/dens**2

@add_quantity(name='gamc_s', args='DT', dependencies=['P{s}_DT', 'P{s}_DU{s}'])
def isentropic_gamc_s(tab, spec, *XYf):
    """Isentropic sound speed"""
    XYf_DT = list(XYf)
    dens = XYf[0]
    pres = tab.get_table('P{s}_DT', spec)(*XYf_DT)

    return dens*tab.q['Cs2', spec](*XYf_DT)/pres

@add_quantity(name='Cs2_check', args='DT', dependencies=['St_DPt', 'Pt_DT'], spec=['t'])
def isentropic_sound_speed2_check(tab, spec, *XYf):
    """Isentropic sound speed square: alternative method."""
    XYf = list(XYf)
    XYf[1] = tab.Pt_DT(*XYf) # now the independant variables are DPt
    return -tab.St_DPt.dFx(*XYf)/tab.St_DPt.dFy(*XYf)

@add_quantity(name='beta_s', args='DT', dependencies=['Cs2'])
def isentropic_beta(tab, spec, *XYf):
    """Isentropic bulk modulus"""
    return XYf[0]*tab.q['Cs2', spec](*XYf)

@add_quantity(name='Ks', args='DT', dependencies=['beta_s'])
def isentropic_compressibility(tab, spec, *XYf):
    """Isentropic compressibility"""
    return 1./(tab.q['beta_s', spec](*XYf))



# Isothermal quantities
# ---------------------
@add_quantity(name='Ct2', args='DT', dependencies=['P{s}_DT'])
def isothermal_sound_speed2(tab, spec, *XYf):
    """Isothermal sound speed square"""
    return tab.get_table('P{s}_DT', spec).dFx(*XYf)

@add_quantity(name='beta_t', args='DT', dependencies=['Ct2'])
def isothermal_beta(tab, spec, *XYf):
    """Isothermal bulk modulus"""
    return XYf[0]*tab.q['Ct2', spec](*XYf)

@add_quantity(name='Kt', args='DT', dependencies=['beta_t'])
def isothermal_compressibility(tab, spec, *XYf):
    """Isothermal compressibility"""
    return 1./(tab.q['beta_t', spec](*XYf))


# Other quantities
# ----------------

@add_quantity(name='gruneisen_coef', args='DT', dependencies=['P{s}_DU{s}', 'U{s}_DT'])
def gruneisen_coefficient(tab, spec, *XYf):
    """Gruneisen Coefficient"""
    XYf = list(XYf)
    XYf[1] = tab.get_table('U{s}_DT', spec)(*XYf)
    return tab.get_table('P{s}_DU{s}', spec)(*XYf)


@add_quantity(name='Cv', args='DT', dependencies=['U{s}_DT'])
def specific_heat_v(tab, spec, *XYf):
    """Constant volume specific heat"""
    return tab.get_table('U{s}_DT',spec).dFy(*XYf)

@add_quantity(name='Cp', args='DT', dependencies=['Cs2', 'Cv', 'Ct2'])
def specific_heat_p(tab, spec, *XYf):
    """Constant pressure specific heat"""
    return tab.q['Cs2', spec](*XYf)*tab.q['Cv', spec](*XYf)\
            /tab.q['Ct2', spec](*XYf)


@add_quantity(name='alpha_exp', args='DT', dependencies=['P{s}_DT', 'Kt'])
def thermal_expansion_coeff(tab, spec, *XYf):
    """Thermal expansion coefficient"""
    return tab.q['Kt', spec](*XYf)*\
            tab.get_table('P{s}_DT', spec).dFy(*XYf)

@add_quantity(name='alpha_exp_check', args='DT', dependencies=['Pt_DT','D_PtT'],
        spec=['t'])
def thermal_expansion_coeff_check(tab, spec, *XYf):
    r"""Thermal expansion coefficient: alternative method
       $\alpha_exp = \frac{1}{V} \left.\frac{\partial V}{\partial T}\right|_P$
    """
    XYf = list(XYf)
    rho = XYf[0].copy()
    XYf[0] = tab.Pt_DT(*XYf)
    return -tab.D_PtT.dFy(*XYf)/rho

#@add_quantity(name='gamc_s', args='DT', dependencies=['P{s}_DT', 'S{s}_DT', 'P{s}_DS{s}'])
#def gamma_adiabatic(tab, spec, *XYf):
#    r"""Adiabatic exponent:
#        $\gamma_S = \frac{\rho}{P}\left.\frac{\partial P}{\partial \rho}\right|_S$"""
#    XYf_DT = list(XYf)
#    rho = XYf[0]
#    P =  tab.get_table('P{s}_DT', spec)(*XYf_DT)
#    XYf_DS = XYf_DT[:]
#    XYf_DS[1] = tab.get_table('S{s}_DT', spec)(*XYf_DT)
#    return rho*tab.get_table('P{s}_DS{s}', spec).dFx(*XYf_DS)/P

@add_quantity(name='gamc_t', args='DT', dependencies=['P{s}_DT'])
def gamma_isothermal(tab, spec, *XYf):
    r"""Isothermal exponent:
        $\gamma_T = \frac{\rho}{P}\left.\frac{\partial P}{\partial \rho}\right|_T$
    """
    XYf_DT = list(XYf)
    rho = XYf[0]
    P =  tab.get_table('P{s}_DT', spec)(*XYf_DT)
    return rho*tab.get_table('P{s}_DT', spec).dFx(*XYf_DT)/P

@add_quantity(name='gamc_s_bad', args='DT', dependencies=['P{s}_DT', 'Ks'])
def gamma_adiabatic_2(tab, spec, *XYf):
    """Adiabatic exponent"""
    return 1./(tab.get_table('P{s}_DT', spec)(*XYf)*tab.q['Ks', spec](*XYf))

@add_quantity(name='Gamma', args='DT', dependencies=['alpha_exp', 'Cv', 'Kt'])
def gruneisen_coefficient_v2(tab, spec, *XYf):
    """Gruneisen Coefficient"""
    return tab.q['alpha_exp'](*XYf)/(XYf[0]*tab.q['Cv'](*XYf)*tab.q['Kt'](*XYf))

@add_quantity(name='g', args='DT', dependencies=['P{s}_DT','Cv'])
def dimensionless_specific_heat(tab, spec, *XYf):
    r"""Dimensionless specific heat:
        $\textrm{g} = \frac{P V}{C_V T}$"""
    return tab.get_table('P{s}_DT', spec)(*XYf)/(XYf[0]*tab.q['Cv'](*XYf)*XYf[1])

@add_quantity(name='G3', args='DT', dependencies=['P{s}_DS{s}', 'gamc_s', 'P{s}_DT', 'S{s}_DT'])
def fundemental_derivative(tab, spec, *XYf):
    r"""Fundemental derivative:
        $\textrm{C} = 1 + \frac{1}{2}\frac{1}{\gamma_s P V^2}\left.\frac{\partial^2 P}{\partial \rho^2}\right|_S$"""
    XYf_DT = XYf[:]
    rho = XYf_DT[0].copy()
    XYf_DS = list(XYf[:])
    XYf_DS[1] = tab.get_table('S{s}_DT', spec)(*XYf_DT)
    return 1 + 0.5*rho**2*tab.get_table('P{s}_DS{s}', spec).dFxx(*XYf_DS)/\
            (tab.q['gamc_s'](*XYf)*tab.get_table('P{s}_DT', spec)(*XYf_DT))

@add_quantity(name='therm_consistency_T', args='DT', dependencies=['U{s}_DS{s}', 'P{s}_DT', 'S{s}_DT'])
def check_thermodynamic_consistency(tab, spec, *XYf):
    """Thermodynamic consistency: relative error"""# ε = (∂E/∂S|_V - T)/T
    XYf_DT = XYf[:]
    T = XYf_DT[1]
    XYf_DS = list(XYf[:])
    XYf_DS[1] = tab.get_table('S{s}_DT', spec)(*XYf_DT)
    # (∂E/∂S|_V - T)/T
    dE_S = tab.get_table('U{s}_DS{s}', spec).dFy(*XYf_DS)
    return np.abs((dE_S - T))/T

@add_quantity(name='therm_consistency_S', args='DT', dependencies=['U{s}_DS{s}', 'P{s}_DT', 'S{s}_DT'])
def check_thermodynamic_consistency(tab, spec, *XYf):
    """Thermodynamic consistency: relative error"""# ε = (∂E/∂V|_S - P)/P )"""
    XYf_DT = XYf[:]
    rho = XYf_DT[0]
    P = tab.get_table('P{s}_DT', spec)(*XYf_DT)
    XYf_DS = list(XYf[:])
    XYf_DS[1] = tab.get_table('S{s}_DT', spec)(*XYf_DT)
    # (∂E/∂V|_S - P)/P
    dE_V = rho**2*tab.get_table('U{s}_DS{s}', spec).dFx(*XYf_DS)
    return  np.abs(dE_V - P)/np.abs(P)

@add_quantity(name='game0', args='DT', dependencies=['P{s}_DT', 'U{s}_DT'])
def game0(tab, spec, *XYf):
    r"""Polytropic index: $P= (\gamma_e - 1) E \rho $"""
    return tab.get_table('P{s}_DT', spec)(*XYf)/(XYf[0]*tab.get_table('U{s}_DT', spec)(*XYf)) + 1

@add_quantity(name='game1', args='DT', dependencies=['P{s}_DU{s}', 'U{s}_DT'])
def game1(tab, spec, *XYf):
    r"""Polytropic index: $P= (\gamma_e - 1) E \rho $"""
    XYf = list(XYf)
    dens = XYf[0]
    XYf[1] = tab.get_table('U{s}_DT', spec)(*XYf)
    return tab.get_table('P{s}_DU{s}', spec).dFy(*XYf)/dens + 1
