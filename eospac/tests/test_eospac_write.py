#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from eospac import EosMaterial
from eospac.eospac.libseswrite import _write_sesbin

import numpy as np

from numpy.testing import assert_allclose

def is_equal(x,y):
    assert x == y


def test_write():
    matid_orid = 3720
    matid_new = 93720
    filename = '/tmp/sesame_test_custom_3720'
    tab_regexp = ['.*t_DT']

    # Reading the original file
    m_ref = EosMaterial(matid_orid, tables=tab_regexp,
            backend='eospac')

    # Writing a custom binary SESAME file
    m_ref.save('/tmp/sesame_test_custom_3720', ext='sesame_bin', matid=matid_new)

    # reading the custom binary file

    m_new = EosMaterial(matid_new, tables=tab_regexp,
                    backend='eospac')


    #print m_ref.Pt_DT['F_Array'].shape
    #print m_new.Pt_DT['F_Array'].shape

    #plt.imshow(np.log(m_new.Pt_DT['F_Array']))
    #plt.show()
    for name in  m_ref.tables:
        if 'S' in name:
            continue
        t_ref, t_new = [o.get_table(name) for o in [m_ref, m_new]]

        for prop_key in ['Mean_Atomic_Num', 'Mean_Atomic_Mass', 'Normal_Density',
                         'R_Array', 'T_Array', 'F_Array']:
            yield assert_allclose, t_ref[prop_key], t_new[prop_key]

