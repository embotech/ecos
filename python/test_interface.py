# for python 3 testing compatibility
from __future__ import print_function
import platform

def import_error(msg):
  print()
  print("## IMPORT ERROR:", msg)
  print()

try:
  from nose.tools import assert_raises, assert_almost_equals
except ImportError:
  import_error("Please install nose to run tests.")
  raise

try:
  import ecos
except ImportError:
  import_error("You must install the ecos module before running tests.")
  raise

try:
  import numpy as np
except ImportError:
  import_error("Please install numpy.")
  raise

try:
  import scipy.sparse as sp
except ImportError:
  import_error("Please install scipy.")
  raise

# global data structures for problem
c = np.array([-1.])
h = np.array([4., -0.])
G = sp.csc_matrix([1., -1.]).T.tocsc()
A = sp.csc_matrix([1.])
b = np.array([3.])
dims = {'q': [], 'l': 2}

def check_solution(solution, expected):
  assert_almost_equals(solution, expected, places=5)

def test_problems():
  myopts = {'feastol': 2e-8, 'reltol': 2e-8, 'abstol': 2e-8, 'verbose':True};
  sol = ecos.solve(c, G, h, dims, **myopts)
  yield check_solution, sol['x'][0], 4

  sol = ecos.solve(c, G, h, dims, A, b, **myopts)
  yield check_solution, sol['x'][0], 3

  new_dims = {'q':[2], 'l': 0}
  sol = ecos.solve(c, G, h, new_dims, **myopts)
  yield check_solution, sol['x'][0], 2

if platform.python_version_tuple() < ('3','0','0'):
  def test_problems_with_longs():
    new_dims = {'q': [], 'l': long(2)}
    myopts = {'feastol': 2e-8, 'reltol': 2e-8, 'abstol': 2e-8};
    sol = ecos.solve(c, G, h, new_dims, **myopts)
    yield check_solution, sol['x'][0], 4

    sol = ecos.solve(c, G, h, new_dims, A, b, **myopts)
    yield check_solution, sol['x'][0], 3

    new_dims = {'q':[long(2)], 'l': 0}
    sol = ecos.solve(c, G, h, new_dims, **myopts)
    yield check_solution, sol['x'][0], 2

def check_keyword(error_type, keyword, value):
  assert_raises(error_type, ecos.solve, c,G,h,dims, **{keyword: value})

def test_failures():
  yield assert_raises, TypeError, ecos.solve
  yield assert_raises, TypeError, ecos.solve, c, G, h, dims, A

  yield assert_raises, ValueError, ecos.solve, c, G, h, {'q':[], 'l':0}
  yield assert_raises, TypeError, ecos.solve, c, G, h, {'q':[4], 'l':-2}

  yield check_keyword, TypeError, 'verbose', 0
  yield check_keyword, ValueError, 'feastol', 0
  yield check_keyword, ValueError, 'abstol', 0
  yield check_keyword, ValueError, 'reltol', 0
  yield check_keyword, ValueError, 'feastol_inacc', 0
  yield check_keyword, ValueError, 'abstol_inacc', 0
  yield check_keyword, ValueError, 'reltol_inacc', 0
  yield check_keyword, ValueError, 'max_iters', -1
  yield check_keyword, TypeError, 'max_iters', 1.1
