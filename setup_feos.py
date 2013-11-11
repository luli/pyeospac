#!/usr/bin/python
# -*- coding: utf-8 -*-

from setuptools import Extension
import os.path
import re
from datetime import datetime
import numpy

def setup_feos(cfg):

    FEOS_INCLUDE = cfg['path']
    FEOS_LIB = cfg['path']

    if not os.path.exists(FEOS_INCLUDE):
        raise OSError("Path does not exist: '{0}'. Please edit setup.cfg !".format(FEOS_INCLUDE))



    #===============================================================================#

    return  [Extension("eospac.feos.libpyfeos",
                 sources=["eospac/feos/libpyfeos.pyx"],
                 include_dirs=[numpy.get_include(), FEOS_INCLUDE],
                 #library_dirs=[FEOS_LIB],
                 #libraries=['feos'],
                 extra_compile_args=['-g'],
                 language='c',
                 extra_link_args=["-lstdc++"],
                 extra_objects=[os.path.join(FEOS_LIB, 'libfeos.a')]
                 )]




