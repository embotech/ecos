import ecos
#import cvxopt as o

from guppy import hpy
# A = o.spmatrix([1.,1.], [0,1],[0,1])
# b = o.matrix([1.,1.])
# G = o.spmatrix([-1.,-1.], [0,1], [0,1])
# h = o.matrix([0.,0.])
# c = o.matrix([1.,1.])
# dims = {'l': 2}

import numpy as np
from scipy import sparse
#ij = np.array([[0,1],[0,1]])

D = np.array([[2.,1.], [3.,4.], [-1.,0.], [0.,-1.]])

#A = sparse.csc_matrix(([1.,1.], ij), (2,2))
#b = np.array([1.,1.])
G = sparse.csc_matrix(D)
h = np.array([4.,12.,0.,0.])
c = np.array([-1.,-1.])
dims = {'l': 4}

sol = ecos.solve(c,G,h,dims) #,A,b)

print sol
print sol['x']

print c
print h
print D
# Ai = numpy.matrix([0,1])
# Ap = numpy.matrix([0,1,2])
# b = numpy.matrix([1.,1.])
# c = numpy.matrix([1.,1.])
# 
# solution = p.solve(Ax, Ai, Ap, b, c)
# print solution['x']
# print solution['y']
# print solution['status']
# 
# print getrefcount(solution['x'])

#h = hpy()
#print h.heap()