import _ecos
from warnings import warn
import numpy as np
from scipy import sparse

def solve(c,G,h,dims,A=None,b=None,verbose=True):
    """ This Python routine "unpacks" scipy sparse matrices G and A into the
        data structures that we need for calling ECOS' csolve routine.
        
        If G and h are both None, then we will automatically create an "empty"
        CSC matrix to use with ECOS.
        
        It is *not* compatible with CVXOPT spmatrix and matrix, although
        it would not be very difficult to make it compatible. We put the
        onus on the user to convert CVXOPT matrix types into numpy, scipy
        array types.
    """
    if G is not None and not sparse.issparse(G):
        raise TypeError("G is required to be a sparse matrix")
    if A is not None and not sparse.issparse(A):
        raise TypeError("A is required to be a sparse matrix")
    
    if G is not None and not sparse.isspmatrix_csc(G):
        warn("Converting G to a CSC matrix; may take a while.")
        G = G.tocsc()
    if A is not None and not sparse.isspmatrix_csc(A):
        warn("Converting A to a CSC matrix; may take a while.")
        A = A.tocsc()
    
    if G is None: m,n1 = 0,len(c)
    else: m,n1 = G.shape
    if A is None: p,n2 = 0,n1
    else: p,n2 = A.shape
    
    if n1 != n2:
        raise TypeError("Columns of A and G don't match")
        
        
    # G.sort_indices() # ECHU: performance hit? do we need this?
    # if A is not None: A.sort_indices()
    
    if G is None:
        if h is not None: raise TypeError("G and h must be supplied together")
        data = np.zeros((0,),dtype=np.double)
        indices = np.zeros((0,),dtype=np.int)
        colptr = np.zeros((n1+1,),dtype=np.int)
        h = np.zeros((0,))
    else:
        data, indices, colptr = G.data, G.indices, G.indptr
    
    
    if A is None:
        if b is not None: raise TypeError("A and b must be supplied together")
        return _ecos.csolve((m,n1,p), c, data, indices, colptr, h, dims, verbose=verbose)
    else:
        return _ecos.csolve((m,n1,p), c, data, indices, colptr, h, dims, A.data, A.indices, A.indptr, b, verbose)
