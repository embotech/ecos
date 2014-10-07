import ecos
import numpy as np
from scipy import *
import scipy.sparse as sp

c = np.array([-1., -1.])
h = np.array([ 4., 12.])
bool_idx = np.array([1])
G = sp.csc_matrix( (array([2.0, 3.0, 1.0, 4.0]), 
	array([0, 1, 0, 1]), 
	array([0, 2, 4])) )  

#idxint n = 2;
#idxint m = 2;
#pfloat feas_Gx[4] = {2.0, 3.0, 1.0, 4.0};
#idxint feas_Gp[3] = {0, 2, 4};
#idxint feas_Gi[4] = {0, 1, 0, 1};

#pfloat feas_c[2] = {-1., -1.};
#pfloat feas_h[2] = {4., 12.};

#idxint bool_idx[1] = {1};

dims = dict()
dims['l'] = 2

sol = ecos.solve(c, G, h, dims, bool_idx=bool_idx)

print sol['x']