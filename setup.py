#!/usr/bin/python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
import os.path
from Cython.Distutils import build_ext
from ConfigParser import ConfigParser
from setup_eospac import setup_eospac

cfg_obj = ConfigParser()
cfg_obj.read('setup.cfg')
cfg_eospac = dict(cfg_obj.items('DEFAULT'))
cfg_feos = dict(cfg_obj.items('FEOS'))

warning_sign =   ''.join([ "="*80, "\n",
                           ' '*10, '{0}',
                           "="*80, "\n"])
extensions_found = []
if os.path.exists(cfg_eospac['path']):
    extensions_found.append(setup_eospac(cfg_eospac))
else:
    print warning_sign.format('Warning: linking to the EOSPAC library was not done. ')





setup(
    name="eospac",
    version='0.1',
    author='Roman Yurchak (LULI)',
    author_email='rth@crans.org',
    packages=find_packages(),
    cmdclass={'build_ext': build_ext},
    ext_modules=[
                 ],
    entry_points = {
        'console_scripts': [ 'eostab2vtk = eospac.viz:convert_to_vtk' ] },
    tests_require=['nose']
)


