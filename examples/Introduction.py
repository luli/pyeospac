# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# PyEospac: exemples of use
# ========================

# <codecell>

import eospac as eos
import numpy as np

# <markdowncell>

# I. Avalable Tables, Information and Options
# --------------------------------------------
# are automatically generated from `eos_Interface.h` at install time

# <codecell>

print eos.tables.keys()

# <codecell>

print eos.info.keys()

# <codecell>

print eos.options.keys()

# <markdowncell>

# II. Single material interpolations
# ----------------------------------

# <codecell>

# Load some tables for the material no 3720 (Aluminium)
# options keys accepts regular expressions. For instance
#    - '.*' would match any table
#    - '[PU]t_DT' matches 'Pt_DT', 'Ut_DT'
#    - etc..
mat = eos.EosMaterial(3720, ['Pt_DT', 'Ut_DT','St_DT', 'Pt_DSt'],
                      options={'.*': {'create_tzero': True},
                               '[PU]t_DT': {'x_convert': 1.0},
                               'St_DT': {'linear': True},
                                })

# see what options are applied
for table_name in mat.tables:
    print table_name, mat.get_table(table_name).options

# <codecell>

# Separate tables can be accessed as parameters of the 'mat' object for instance 'm.Pt_DT'
# this is a shortcut from m.get_table('Pt_DT')
# Info values can be acessed as attributes of the EosMaterial class
print mat.Pt_DT['Normal_Density'], mat.Pt_DT['Mean_Atomic_Num']

# <codecell>

# The list of all avalables info parameters for a given table can be obtained with
print mat.Pt_DT

# <codecell>

X = np.array([0.1, 0.2, 0.3]) # density array in g.cm⁻³
Y = np.array([1e5, 1e6, 1e8]) # temperature array in K
# interpolates the total pressure at given X,Y points
# X,Y could be arrrays of any dimension
print 'F', mat.Pt_DT(X,Y)

# ∂F/∂x  first order derivative as returned by EOSPAC6 (same thing for dFy)
print 'dFx', mat.Pt_DT.dFx(X,Y)

# ∂²F/∂x² second order derivative (finite differences over EOSPAC's first derivative)
print 'dFxx', mat.Pt_DT.dFxx(X,Y)

# ∂²F/∂x∂y second order derivative 
print 'dFxy', mat.Pt_DT.dFxy(X,Y)
print 'dFyx', mat.Pt_DT.dFyx(X,Y)
# checking that F is C² : ||dFxy - dFyx|| << 1
print 'max(dFxy - dFxy)', np.abs(mat.Pt_DT.dFxy(X,Y) - mat.Pt_DT.dFyx(X,Y)).max()

# <markdowncell>

# III. Mixtures
# -------------

# <codecell>

# The EosMixture class allows to manipulate multiple matierial with a similar interface to EosMaterial's one
# The following loads tables Pt_DT, Ut_DT, St_DT and Pt_DSt for materials 3720 and 7593 using regular expressions.
# The setup options are very similar to EosMaterial class. An additional possibility is to provide a tuple as key
# in order to specify the material the option should be applied to, otherwise the corresponding option will be applyed
# to all materials
mix = eos.EosMixture([3720, 7593], ['[PUS]t_DT', 'Pt_DSt'],
                      options={'.*': {'create_tzero': True},
                               '[PU]t_DT': {'x_convert': 1.0},
                               (7593, 'St_DT'): {'linear': False},
                                })
print 'Loaded tables:', mix.tables

# <codecell>

# It is possible access to different loaded materials separately 
print "object:", mix.m3720

# See previous section to all possibilities of EosMaterial class. For instance:
print 'interpolating :', mix.m3720.Pt_DT(X,Y) 

# <codecell>

# Multi-material properties can be accessed by specifying the table name as attribute. For instance:
print mix.Pt_DT['Material_ID']
print mix.Pt_DT['Normal_Density']

# <codecell>

# Assuming that X, Y used in interpolations are of shape (N1,N2,N3, etc..)
# the shape of the fraction array used by eos_Mix should be (number_of_materials, N1, N2, N3, etc..)
# Let's make a 25/75 mixture of materials 7320 and 7593
frac = np.zeros((2,) + X.shape)
frac[0] = 0.25
frac[1] = 1 - frac[0]
print 'frac.shape', frac.shape
print 'X.shape', X.shape

# <codecell>

# the interpolation within mixture is then done as follows:

print 'F', mix.Pt_DT(X,Y, frac)
print 'dFx', mix.Pt_DT.dFx(X,Y, frac)

# <markdowncell>

# IV. Derived quantities
# -----------------------
# A mecanism to access derived quanties is provided (inspired from `yt toolkit`). Most of the quantities described in EOSPAC's manual and a few additional ones are implemented.

# <codecell>

# see eospac/quantities.py for more details
print eos.derived_quantities

# <codecell>

# If a derived quantity is given to EosMaterial (or EosMixture), all the required tables to 
# compute the corresponding quantity will be pulled automatically
mat = eos.EosMaterial(3720, ['Pt_DT', 'gamc_s','gamc_t', 'Cp'], spec=['t', 'ic'])
print 'mat.tables', mat.tables
print 'Cp dependencies:', eos.derived_quantities['Cp']['dependencies']
print 'Cp required tables:', eos.derived_quantities.get_dependencies('Cp', spec=['t'])

# <codecell>

# Now we can acess for instance the gamc_t quantity in a similar way to a table
print 'gamc_t total', mat.q['gamc_t','t'](X,Y)
print 'gamc_t ion+cc', mat.q['gamc_t','ic'](X,Y)

# The derived quantities work with mixtures as well:
print 'gamc_t mixture total', mix.q['gamc_t','t'](X,Y, frac)

# <codecell>

# It is very easy to define new derived quantities
# Let's define a quantity gamc_t_v0 which actually does the same thing as gamc_t

@eos.add_quantity(name='gamc_t_v0', args='DT', dependencies=['Pt_DT'])
def gamma_isothermal_v0(tab, spec, rho, temp):
    r"""Isothermal gamma exponent"""
    P =  tab.Pt_DT(rho, temp)
    return rho*tab.Pt_DT.dFx(rho, temp)/P

print 'gamc_t_v0 total', mat.q['gamc_t_v0','t'](X,Y)

# <codecell>

# The previous definition is rather intuitive, however it is restrictive as it will only 
# work for total tables and a single material.
# A more general definition allow to take into account different species (total, electron, ion+cc, etc)
# as well as multi-material setups.
# In this case 'P{s}_DT' will be expanded to 'Pt_DT','Pe_DT', etc, according to the value of spec.
# spec could be 't', 'e', 'ie', etc..
# *XYf is a tuple containing:
#   - (X,Y) if applied to EosMaterial
#   - (X, Y, frac) if applied to EosMixture


@eos.add_quantity(name='gamc_t_v1', args='DT', dependencies=['P{s}_DT'])
def gamma_isothermal_v1(tab, spec, *XYf):
    r"""Isothermal gamma exponent """
    rho = XYf[0]
    P =  tab.get_table('P{s}_DT', spec)(*XYf)
    return rho*tab.get_table('P{s}_DT', spec).dFx(*XYf)/P

print 'gamc_t_v1 total', mat.q['gamc_t_v1','t'](X,Y)
print 'gamc_t_v1 mixture total', mix.q['gamc_t_v1','t'](X,Y, frac)

# <markdowncell>

# V. Data visualization
# ---------------------
# a few examples of how to plot EoS data with `matplotlib`

# <codecell>

rc= plt.rcParams
rc['figure.figsize'] = (12, 8)
rc['figure.dpi'] = 300
rc['font.size'] = 20

# <codecell>

mat = eos.EosMaterial(3720, ['Pt_DT', 'Ut_DT','gamc_s', 'therm_consistency'], spec=['t'])

# <codecell>

eos.eos_plot(mat, 'Pt_DT', ax=plt.subplot(111))

# <codecell>

# the white contour corresponds to γ_T = 1
eos.eos_plot(mat, 'gamc_t', ax=plt.subplot(111))

# <codecell>

# Checking that dE =  TdS - PdV
eos.eos_plot(mat, 'therm_consistency', ax=plt.subplot(111), vmax=1.0)

