from distutils.core import setup, Extension
from glob import glob

ecos = Extension('ecos', #libraries = ['m'],
                    # define LDL and AMD to use long ints
                    define_macros = [('DLONG', None), ('LDL_LONG', None)],
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