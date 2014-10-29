import ecos
import numpy as np
from scipy import *
import scipy.sparse as sp

c = np.array([-1., -1.])
h = np.array([ 4., 12.])
bool_idx = np.array([float(1.0) ])
G = sp.csc_matrix( (array([2.0, 3.0, 1.0, 4.0]), 
	array([0, 1, 0, 1]), 
	array([0, 2, 4])) )  

dims = dict()
dims['l'] = 2

sol = ecos.solve(c, G, h, dims, integer_vars_idx=bool_idx)

print sol['x']

c = np.array([-1., -1.])
h = np.array([ 4., 12.])
bool_idx = np.array([float(0.0)])
G = sp.csc_matrix( (array([2.0, 3.0, 1.0, 4.0]), 
	array([0, 1, 0, 1]), 
	array([0, 2, 4])) )  

dims = dict()
dims['l'] = 2

sol = ecos.solve(c, G, h, dims, integer_vars_idx=bool_idx)

print sol['x']