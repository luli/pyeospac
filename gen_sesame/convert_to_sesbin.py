#!/usr/bin/python
# -*- coding: utf-8 -*-

import os.path
import eospac as eos
from eospac.tabulated import avalable_tabs

NESTED_DIR = "/home/rth/luli/NestedOutflows/NestedOutflows/simulations/NestedOutflows"
DATA_DIR = "/home/rth/src/pyeospac/data/"



#==============================================================================
#                             Aluminum
#==============================================================================
#mat1 = eos.EosMaterial(73719,
#      options={'type': 'ionmix',
#          'abar':  26.9815, 'zbar':13.0, 'rho_ref': 2.7,
#          'create_tzero': 'linear',
#              'path': os.path.join(NESTED_DIR, 'al-imx-32g.cn4')},
#      units='cgs', backend='tabulated')
#mat1.save(os.path.join(DATA_DIR, 'Al_ionmix_73719.sesb'))
#
#feos_opts ={'use_maxwell': True, 'use_softspheres': True,
#        'maxwell_temp_arr': None, 'max_materials': 2,
#        'grid_subsample_default': 0}
#
#mat2 = eos.EosMaterial(3717, tables=[tab for tab in avalable_tabs if 't' in tab],
#      options=feos_opts,
#      units='cgs', backend='feos')
#mat2.save(os.path.join(DATA_DIR, 'Al_feos_83719.sesb'), matid=83719)

#==============================================================================
#                             Polystyrene
#==============================================================================
#matbase = 7592
#material = 'polystyrene'
#tab_options={'type': 'ionmix',
#          'abar':  6.51, 'zbar':3.5, 'rho_ref': 1.044,
#          'create_tzero': 'linear',
#          'path': 'polystyrene-imx-32g.cn4'}
#
#feos_opts ={'use_maxwell': True, 'use_softspheres': False,
#        'maxwell_temp_arr': None, 'max_materials': 2,
#        'grid_subsample_default': 0}
#==============================================================================
#                             Aluminum
#==============================================================================
matbase = 3719
feosid = 83719
material = 'Al'
#tab_options={'type': 'ionmix',
#          'abar':  26.9815, 'zbar':13, 'rho_ref': 2.7,
#          'create_tzero': 'linear',
#          'path': 'al-imx-32g.cn4'}

feos_opts ={'use_maxwell': True, 'use_softspheres': True,
        'maxwell_temp_arr': None, 'max_materials': 1,
        'grid_subsample_default': 0}
#==============================================================================
#                               Iron
#==============================================================================
#matbase = 2140
#material = 'Fe'
#tab_options={'type': 'ionmix',
#          'abar':  55.845, 'zbar':26, 'rho_ref': 7.85,
#          'create_tzero': 'linear',
#          'path': 'fe-imx-32g.cn4'}
#
#feos_opts ={'use_maxwell': True, 'use_softspheres': True,
#        'maxwell_temp_arr': None, 'max_materials': 2,
#        'grid_subsample_default': 0}

#tab_options['path'] = os.path.join(NESTED_DIR, tab_options['path'])
#mat1 = eos.EosMaterial(7e4+matbase,
#      options=tab_options,
#      units='cgs', backend='tabulated')
#mat1.save(os.path.join(DATA_DIR, '{0}_ionmix_{1}.sesb').format(material, int(7e4+matbase)))

#==============================================================================
#                               Iron
#==============================================================================
#matbase = 3332
#feosid = 83333
#material = 'Cu'
#
#feos_opts ={'use_maxwell': False, 'use_softspheres': True,
#        'maxwell_temp_arr': None, 'max_materials': 1,
#        'grid_subsample_default': 1}
#
mat2 = eos.EosMaterial(matbase, tables=[tab for tab in avalable_tabs],
      options=feos_opts,
      units='cgs', backend='feos')
if feos_opts['use_maxwell']:
    print 'Critical point {3}: {0:.3f} g.cm⁻³,  {1:.1f} K,  {2:.0f} bar'.format(
        mat2.Pt_DT['critical_rho'], mat2.Pt_DT['critical_temp']*11640,
        mat2.Pt_DT['critical_P']*1e-10*1e4, material)
mat2.save(os.path.join(DATA_DIR, '{0}_feos_{1}.sesb'.format(material, feosid)), matid=feosid)
