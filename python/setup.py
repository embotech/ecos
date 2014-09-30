from __future__ import print_function
try:
    from setuptools import setup, Extension
except ImportError:
    print("Please use pip (https://pypi.python.org/pypi/pip) to install.")
    raise

from glob import glob
from platform import system
import numpy

lib = []
if system() == 'Linux':
    lib += ['rt']

_ecos = Extension('_ecos', libraries = lib,
                    # define LDL and AMD to use long ints
                    # also define that we are building a python module
                    define_macros = [
                        ('PYTHON',None),
                        ('DLONG', None),
                        ('LDL_LONG', None)],
                    include_dirs = ['../include', numpy.get_include(),
                        '../external/amd/include',
                        '../external/ldl/include',
                        '../external/SuiteSparse_config'],
                    sources = ['ecosmodule.c',
                        '../external/ldl/src/ldl.c'
                    ] + glob('../external/amd/src/*.c')
                      + glob('../src/*.c'))

setup(
    name = 'ecos',
    version = '1.0.5',
    author = 'Alexander Domahidi, Eric Chu',
    author_email = 'domahidi@embotech.com, echu@cs.stanford.edu',
    url = 'http://github.com/ifa-ethz/ecos',
    description = 'This is the Python package for ECOS: Embedded Cone Solver. See Github page for more information.',
    license = "GPLv3",
    py_modules = ['ecos'],
    ext_modules = [_ecos],
    install_requires = [
        "numpy >= 1.7",
        "scipy >= 0.12"
    ]
)
