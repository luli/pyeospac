#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from eospac import EosMaterial
from eospac.eospac.libsesio import _write_sesbin

import numpy as np

from numpy.testing import assert_allclose

def is_equal(x,y):
    assert x == y


