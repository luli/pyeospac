#PyEOSPAC: Python bindings for the EOSPAC6 library

## Overview



## Dependencies

Python 2.7, `setuptools`, `cython`, `numpy`, `scipy` modules are necessary for the base functionallity.

Besides, at least one of the following libraries should also be installed,

 * the [EOSPAC library](https://laws.lanl.gov/projects/data/eos/eospacReleases.php) (Pimentel et al., LANL)
 * the [FEOS code](http://th.physik.uni-frankfurt.de/~faik/feos.php?lang=eng) (Faik et al.): 

otherwise the build-in gamma law EoS will be used. 

Optional dependencies  include `nose`, `matplotlib`, `ipython` and `pandas`.

## Installation notes

1. Edit `eospac6.cfg` file and set the correct the path to EOSPAC6 install.

2. To install without root permissions run


     $ python setup.py install --user


3. Optional: to execute the test suite run,


     $ nosetests ./eospac/

Warning:
 a. this runs some tests shipped with PyEOSPAC, not the EOSPAC6 test suite..
 b. the test suite requires access to tables 3720 and 7593, and therefore a Sesame data files (not included).


## Examples of use

See the `examples/` directory for some examples of use, and in particular https://luli.github.io/pyeospac/Introduction.html
