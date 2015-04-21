#PyEOSPAC: Equations of State (EoS) in Python
## Overview

Initially developed as python wrapper around the EOSPAC library, PyEOSPAC currently supports several backends (EOSPAC, FEOS code) and analytical EoS models (ideal gas, Van der Vaals gas), with an aim to provide a unified and user-friendly access to the underlying Equation of State (EoS) data.

The list of features include,

  * High level access to the EOSPAC routines for tabulated EoS, including interpolation, material mixing and table inversion. 
  * Runtime definition of derived quantities (sound speed, heat capacity, etc.) from higher order derivatives.
  * Calculation of Hugoniot and Isentropes curves 
  * Maxwell constructions. Calculation of binodal and spinodal curves for tables out of equilibrium (i.e. with Van der Waals loops).

Use cases,

  * General EoS visualisation
  * Coupled with the [opacplot2](https://github.com/rth/opacplot2), this module can be used to resample EoS tables, on a user defined grid. Thermodynamic consistency and interpolation error can also be estimated.
  * In the context of hydrodynamic codes, this module can be loaded by the visualization software (e.g. VisIt, `yt`) for visualization or post-processing purposes. 

## Dependencies

Python 2.7, `setuptools`, `cython`, `numpy`, `scipy` modules are necessary for the basic functionality.

Besides, at least one of the following libraries should also be installed,

 * the [EOSPAC library](https://laws.lanl.gov/projects/data/eos/eospacReleases.php)
 * the [FEOS code](http://th.physik.uni-frankfurt.de/~faik/feos.php?lang=eng)

otherwise the build-in analytical EoS models will be used. 

Optional dependencies  include `nose`, `matplotlib`, `ipython` and `pandas`.

## Installation notes

1. Edit `setup.cfg` file and set the correct the path to the EOSPAC6 and/or FEOS install.

2. To install without root permissions run


         $ python setup.py install --user


3. Optional: the test suite can be executed with,


         $ nosetests ./eospac/

Warnings:

  a. this runs the tests shipped with PyEOSPAC, not the EOSPAC6 test suite.
  b. the test suite requires access to tables 3720 and 7593, and therefore a Sesame data files (not included).


## Examples

See the `examples/` directory for some examples of use, and in particular  [Introduction.html](https://luli.github.io/pyeospac/Introduction.html)

## Future work

  * add automatic numerical differentiation with `numdifftools` when calculating the partial derivatives of the EoS.
