from distutils.core import setup, Extension
from glob import glob
from numpy import get_include

direct = Extension('pdos_direct', libraries = ['m'],
                    include_dirs = ['../', get_include(),
                        '../direct/external/AMD/Include', 
                        '../direct/external/LDL/Include',
                        '../direct/external/SuiteSparse_config'],
                    sources = ['pdosmodule.c',
                        '../direct/private.c',
                        '../direct/external/LDL/Source/ldl.c',
                        '../cones.c', '../cs.c', '../pdos.c', '../util.c'
                    ] + glob('../direct/external/AMD/Source/*.c'))

indirect = Extension('pdos_indirect', libraries = ['m'],
                    include_dirs = ['../', get_include()],
                    define_macros = [('INDIRECT', None)],
                    sources = ['pdosmodule.c',
                        '../indirect/private.c',
                        '../cones.c', '../cs.c', '../pdos.c', '../util.c'
                    ])


setup(  name = 'pdos',
        version = '1.0',
        description = 'This is Python package to wrap our first-order solvers',
        ext_modules = [direct, indirect],
        requires = ["numpy (>= 1.7)"])