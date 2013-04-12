from distutils.core import setup, Extension
from glob import glob
from platform import system

lib = []
if system() == 'Linux':
    lib += ['rt']

ecos = Extension('ecos', libraries = lib,
                    # define LDL and AMD to use long ints
                    # also define that we are building a python module
                    define_macros = [
                        ('PYTHON',None),
                        ('DLONG', None),
                        ('LDL_LONG', None)],
                    include_dirs = ['../include',
                        '../external/amd/include', 
                        '../external/ldl/include',
                        '../external/SuiteSparse_config'],
                    sources = ['ecosmodule.c',
                        '../external/ldl/src/ldl.c'
                    ] + glob('../external/amd/src/*.c')
                      + glob('../src/*.c'))


setup(  name = 'ecos',
        version = '1.0',
        description = 'This is Python package for ECOS: Embedded Cone Solver.',
        ext_modules = [ecos],
        requires = ["cvxopt (>= 1.1.5)"])
