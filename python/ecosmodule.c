#include <Python.h>
#include "ecos.h"
#include "cvxopt.h"

static PyObject *ecos(PyObject* self, PyObject *args)
{
  /* Expects a function call 
   *     sol = ecos(c,G,h,dims,A,b)
   * where
   *
   * `c` is a cvxopt (dense) column vector
   * `G` is a cvxopt (sparse) matrix
   * `h` is a cvxopt (dense) column vector
   * `dims` is a dictionary with
   *    `dims['l']` an integer specifying the dimension of positive orthant cone
   *    `dims['q']` an *array* specifying dimensions of second-order cones
   * `A` is an optional argument, which is a cvxopt (sparse) matrix
   * `b` is an optional argument, which is a cvxopt (dense) column vector
   *
   * This call will solve the problem
   *
   *    minimize     c'*x
   *    subject to   A*x = b
   *                 h - G*x \in K
   *
   * The code returns a Python dictionary with five keys, 'x', 'y', 'info', 's',
   * and 'z'. These correspond to the following:
   * 
   * `x`: primal variables
   * `y`: dual variables for equality constraints
   * `s`: slacks for Gx + s <= h, s \in K
   * `z`: dual variables for inequality constraints s \in K
   * `info`: another dictionary with the following fields:
   *    exitflag: 0=OPTIMAL, 1=PRIMAL INFEASIBLE, 2=DUAL INFEASIBLE, -1=MAXIT REACHED
   *  infostring: gives information about the status of solution
   *       pcost: value of primal objective
   *       dcost: value of dual objective
   *        pres: primal residual on inequalities and equalities
   *        dres: dual residual
   *        pinf: primal infeasibility measure
   *        dinf: dual infeasibility measure
   *     pinfres: NaN
   *     dinfres: 3.9666e+15
   *         gap: duality gap
   *      relgap: relative duality gap
   *          r0: ???
   *      numerr: numerical error?
   *        iter: number of iterations
   *      timing: dictionary with timing information
   */

  /* data structures for arguments */
  matrix *c, *h, *b = NULL;
  spmatrix *G, *A = NULL;
  PyObject *dims;
  idxint m,n;
  idxint p = 0;  // number of equality constraints
  
  /* ECOS data structures */
  idxint l = 0;
  idxint *q = NULL;
  idxint ncones = 0;   
  
  pfloat *Gpr = NULL;
  idxint *Gjc = NULL;
  idxint *Gir = NULL;
 
  pfloat *Apr = NULL;    
  idxint *Ajc = NULL;
  idxint *Air = NULL;   
  
  pfloat *cpr = NULL;
  pfloat *hpr = NULL;
  pfloat *bpr = NULL;
  
  pwork* mywork;
  
  idxint i;

  // parse the arguments and ensure they are the correct type 
  if( !PyArg_ParseTuple(args, "OOOO!|OO",
      &c,
      &G,
      &h,
      &PyDict_Type, &dims,
      &A,
      &b)
    ) { return NULL; }
  
  /* set G */
  if ((SpMatrix_Check(G) && SP_ID(G) != DOUBLE)){
      PyErr_SetString(PyExc_TypeError, "G must be a sparse 'd' matrix");
      return NULL;
  }
  if ((m = SP_NROWS(G)) <= 0) {
      PyErr_SetString(PyExc_ValueError, "m must be a positive integer");
      return NULL;
  }
  if ((n = SP_NCOLS(G)) <= 0) {
      PyErr_SetString(PyExc_ValueError, "n must be a positive integer");
      return NULL;
  }
  Gpr = SP_VALD(G);
  Gir = SP_ROW(G);
  Gjc = SP_COL(G);
  
  /* set c */
  if (!Matrix_Check(c) || MAT_NCOLS(c) != 1 || MAT_ID(c) != DOUBLE) {
      PyErr_SetString(PyExc_TypeError, "c must be a dense 'd' matrix with one column");
      return NULL;
  }
  
  if (MAT_NROWS(c) != n){
      PyErr_SetString(PyExc_ValueError, "c has incompatible dimension with G");
      return NULL;
  }
  cpr = MAT_BUFD(c);

  /* set h */
  if (!Matrix_Check(h) || MAT_NCOLS(h) != 1 || MAT_ID(h) != DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "h must be a dense 'd' matrix with one column");
    return NULL;
  }
  
  if (MAT_NROWS(h) != m){
      PyErr_SetString(PyExc_ValueError, "h has incompatible dimension with G");
      return NULL;
  }
  hpr = MAT_BUFD(h);
  
  /* get dims['l'] */
  PyObject *linearObj = PyDict_GetItemString(dims, "l");
  if(linearObj) {
    if(PyInt_Check(linearObj) && ((l = (idxint) PyInt_AsLong(linearObj)) >= 0)) {
      // do nothing
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['l'] ought to be a nonnegative integer");
      return NULL;
    }
  }
  
  /* get dims['q'] */
  PyObject *socObj = PyDict_GetItemString(dims, "q");
  if(socObj) {
    if (PyList_Check(socObj)) {
      ncones = PyList_Size(socObj);
      q = calloc(ncones, sizeof(idxint));
      for (i = 0; i < ncones; ++i) {
          PyObject *qi = PyList_GetItem(socObj, i);
          if(PyInt_Check(qi) && ((q[i] = (idxint) PyInt_AsLong(qi)) > 0)) {
            // do nothing
          } else {
            PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list of positive integers");
            return NULL;
          }

      }
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list");
      return NULL;
    }
  }
  
  if(A && b) {
    /* set A */
    if ((SpMatrix_Check(A) && SP_ID(A) != DOUBLE)){
        PyErr_SetString(PyExc_TypeError, "A must be a sparse 'd' matrix");
        if(q) free(q);
        return NULL;
    } 
    if ((p = SP_NROWS(A)) <= 0) {
        PyErr_SetString(PyExc_ValueError, "m must be a positive integer");
        if(q) free(q);
        return NULL;
    }
    if (SP_NCOLS(A) != n) {
        PyErr_SetString(PyExc_ValueError, "A has incompatible dimension with c");
        if(q) free(q);
        return NULL;
    }
    Apr = SP_VALD(A);
    Air = SP_ROW(A);
    Ajc = SP_COL(A);
     
    /* set b */
    if (!Matrix_Check(b) || MAT_NCOLS(b) != 1 || MAT_ID(b) != DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "b must be a dense 'd' matrix with one column");
        if(q) free(q);
        return NULL;
    }
    if (MAT_NROWS(b) != p){
        PyErr_SetString(PyExc_ValueError, "b has incompatible dimension with A");
        if(q) free(q);
        return NULL;
    }
    bpr = MAT_BUFD(b);
  } else if (A || b) {
    // check that A and b are both supplied
    PyErr_SetString(PyExc_ValueError, "A and b arguments must be supplied together");
    if(q) free(q);
    return NULL;
  }
  
  /* This calls ECOS setup function. */
  mywork = ECOS_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, cpr, hpr, bpr);
  if( mywork == NULL ){
      PyErr_SetString(PyExc_RuntimeError, "Internal problem occurred in ECOS while setting up the problem.\nPlease send a bug report with data to Alexander Domahidi.\nEmail: domahidi@control.ee.ethz.ch");
      if(q) free(q);
      return NULL;
  }
  
  
  /* Solve! */    
  idxint exitcode = ECOS_solve(mywork);
  
  /* create output (all data is *deep copied*) */
  // TODO: request CVXOPT API for constructing from existing pointer
  /* x */
  matrix *x;
  if(!(x = Matrix_New(n,1,DOUBLE)))
    return PyErr_NoMemory();
  for(i = 0; i < n; ++i)
    MAT_BUFD(x)[i] = mywork->x[i];
        
  /* y */
  matrix *y;
  if(!(y = Matrix_New(p,1,DOUBLE)))
    return PyErr_NoMemory();
  for(i = 0; i < p; ++i)
    MAT_BUFD(y)[i] = mywork->y[i];
  
  /* info dict */
  // infostring
  const char* infostring;
  switch( exitcode ){
      case ECOS_OPTIMAL:
          infostring = "Optimal solution found";
          break;
      case ECOS_MAXIT:
          infostring = "Maximum number of iterations reached";
          break;
      case ECOS_PINF:
          infostring = "Primal infeasible";
          break;
      case ECOS_DINF:
          infostring = "Dual infeasible";
          break;
      case ECOS_KKTZERO:
          infostring = "Element of D zero during KKT factorization";
          break;
      case ECOS_OUTCONE:
          infostring = "PROBLEM: Mulitpliers leaving the cone";
          break;
      default:
          infostring = "UNKNOWN PROBLEM IN SOLVER";
  }		
  
  // numerical errors
  idxint numerr = 0;
  if( (exitcode == ECOS_NUMERICS) || (exitcode == ECOS_OUTCONE) || (exitcode == ECOS_FATAL) ){
      numerr = 1;
  }
  
  // timings
#if PROFILING > 0        
	PyObject *tinfos = Py_BuildValue(
#if PROFILING > 1
    "{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
#else
    "{s:d,s:d,s:d}",
#endif
#if PROFILING > 1
    "tkktcreate",(double)mywork->info->tkktcreate,
    "tkktsolve",(double)mywork->info->tkktsolve,
    "tkktfactor",(double)mywork->info->tfactor,
    "torder",(double)mywork->info->torder,
    "ttranspose",(double)mywork->info->ttranspose,
#endif
    "runtime",(double)mywork->info->tsolve + (double)mywork->info->tsetup,
    "tsetup",(double)mywork->info->tsetup,
    "tsolve",(double)mywork->info->tsolve);
#endif

  PyObject *infoDict = Py_BuildValue(
#if PROFILING > 0
    "{s:l,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:l,s:s,s:O,s:l}",
#else
    "{s:l,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:l,s:s,s:l}",
#endif
    "exitFlag", exitcode, 
    "pcost", (double)mywork->info->pcost, 
    "dcost", (double)mywork->info->dcost,
    "pres", (double)mywork->info->pres,
    "dres", (double)mywork->info->dres,
    "pinf", (double)mywork->info->pinf,
    "dinf", (double)mywork->info->dinf,
    "pinfres",(double)mywork->info->pinfres,
    "dinfres",(double)mywork->info->dinfres,
    "gap",(double)mywork->info->gap,
    "relgap",(double)mywork->info->relgap,
    "r0",(double)mywork->stgs->feastol,
    "iter",mywork->info->iter,
    "infostring",infostring,
#if PROFILING > 0
    "timing", tinfos,
#endif
    "numerr",numerr);

  /* s */
  matrix *s;
  if(!(s = Matrix_New(m,1,DOUBLE)))
    return PyErr_NoMemory();
  for(i = 0; i < m; ++i)
    MAT_BUFD(s)[i] = mywork->s[i];
    
  /* z */
  matrix *z;
  if(!(z = Matrix_New(m,1,DOUBLE)))
    return PyErr_NoMemory();
  for(i = 0; i < m; ++i)
    MAT_BUFD(z)[i] = mywork->z[i];
    
  /* cleanup */
  ECOS_cleanup(mywork, 0);
  
  PyObject *returnDict = Py_BuildValue(
    "{s:O,s:O,s:O,s:O,s:O}",
    "x",x,
    "y",y,
    "z",z,
    "s",s,
    "info",infoDict);
    
  return returnDict;
}

static PyMethodDef ECOSMethods[] =
{
  {"ecos", (PyCFunction)ecos, METH_VARARGS, 
    "Solve an SOCP using ECOS."},
  {NULL, NULL, 0, NULL} // sentinel
};

PyMODINIT_FUNC initecos(void)  
{
  PyObject *m;

  m = Py_InitModule("ecos", ECOSMethods);
  
  if (import_cvxopt() < 0) return; // for cvxopt support
  
  if(m == NULL) return;

// #ifdef INDIRECT
//   PDOSError = PyErr_NewException("pdos_indirect.error", NULL, NULL);
// #else
//   PDOSError = PyErr_NewException("pdos_direct.error", NULL, NULL);
// #endif
// 
//   Py_INCREF(PDOSError);
//   PyModule_AddObject(m, "error", PDOSError);
  
//  Py_AtExit(&cleanup);
}
