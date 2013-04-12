import ecos
import cvxopt as o

#from guppy import hpy
A = o.spmatrix([1.,1.], [0,1],[0,1])
b = o.matrix([1.,1.])
G = o.spmatrix([-1.,-1.], [0,1], [0,1])
h = o.matrix([0.,0.])
c = o.matrix([1.,1.])
dims = {'l': 2}

print c
print h

sol = ecos.ecos(c,G,h,dims,A,b)
print sol
print sol['x']
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

# h = hpy()
# print h.heap()