import _ecos
from warnings import warn
import numpy as np
from scipy import sparse

def solve(c,G,h,dims,A=None,b=None):
    """ This Python routine "unpacks" scipy sparse matrices G and A into the
        data structures that we need for calling ECOS' csolve routine.
        
        It is *not* compatible with CVXOPT spmatrix and matrix, although
        it would not be very difficult to make it compatible. We put the
        onus on the user to convert CVXOPT matrix types into numpy, scipy
        array types.
    """
    if not sparse.issparse(G):
        raise "G is required to be a sparse matrix"
    if A is not None and not sparse.issparse(A):
        raise "A is required to be a sparse matrix"
    
    if not sparse.isspmatrix_csr(G):
        warn("Converting G to a CSR matrix; may take a while.")
        G = G.tocsr()
    if A is not None and not sparse.isspmatrix_csr(A):
        warn("Converting A to a CSR matrix; may take a while.")
        A = A.tocsr()
        
    m,n1 = G.shape
    if A is None: p,n2 = 0,n1
    else: p,n2 = A.shape
    
    if n1 != n2:
        raise "Columns of A and G don't match"
        
    # G.sort_indices() # ECHU: performance hit? do we need this?
    # if A is not None: A.sort_indices()
    
    if A is None:
        if b: raise "A and b must be supplied together"
        return _ecos.csolve((m,n1,p), c, G.data, G.indices, G.indptr, h, dims)
    else:
        return _ecos.csolve((m,n1,p), c, G.data, G.indices, G.indptr, h, dims, A.data, A.indices, A.indptr, b)
