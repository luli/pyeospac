PyEOSPAC: Python bindings for EOSPAC6
=====================================

Installation notes
==================

System requirements:

* Python 2.7
* setuptools
* Cython 
* numpy

Note: This module has only been tested on Linux. 

Optional dependencies:

* nose - to run the included test suite
* >=matplotlib-1.0
* ipython

1. Edit `eospac6.cfg` file and set the correct the path to EOSPAC6 install.

2. To install without root permissions run
.. code-block:: bash

    $ python setup.py install --user


3. Optional: to execute the test suite run

.. code-block:: bash

    $ nosetests ./eospac/

Warning: (a) this runs some tests shipped with PyEOSPAC, not the EOSPAC6 test suite.
         (b) the test suite requires access to tables 3720 and 7593, therefore a Sesame data file  


Examples of use
===============

See the `examples/` directory for some examples of use, and in particular `examples/Introduction.html`
