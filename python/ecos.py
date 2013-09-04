import _ecos
import cvxopt

def solve(c,G,h,dims,A=None,b=None):
    c = cvxopt.matrix(c)
    G = cvxopt.sparse(G)
    h = cvxopt.matrix(h)
    if A:
        A = cvxopt.sparse(A)
    if b:
        b = cvxopt.matrix(b)
    return _ecos.csolve(c, G, h, dims, A, b)
