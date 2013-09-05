#include <Python.h>
#include "ecos.h"
#include "numpy/arrayobject.h"

/* IMPORTANT: This code now uses numpy array types. It is a private C module
 * in the sense that end users only see the front-facing Python code in
 * "ecos.py"; hence, we can get away with the inputs being numpy arrays of
 * the CSR data structures.
 *
 * WARNING: This code also does not check that the data for the sparse 
 * matrices are *actually* in column compressed storage for a sparse matrix. 
 * The C module is not designed to be used stand-alone. If the data provided
 * does not correspond to a CSR matrix, this code will just crash inelegantly.
 * Please use the "solve" interface in ecos.py.
 */
//#include "cvxopt.h"

/* ECHU: Note, Python3.x may require special handling for the int and double
 * types. */
static inline int getIntType() {
  switch(sizeof(idxint)) {
    case 1: return NPY_INT8;
    case 2: return NPY_INT16;
    case 4: return NPY_INT32;
    case 8: return NPY_INT64;
    default: return NPY_INT32;  // defaults to 4 byte int
  }
}

static inline int getDoubleType() {
  // ECHU: known bug, if pfloat isn't "double", will cause aliasing in memory
  return NPY_DOUBLE;
}

static inline PyArrayObject *getContiguous(PyArrayObject *array, int typenum) {
  // gets the pointer to the block of contiguous C memory
  // the overhead should be small unless the numpy array has been
  // reordered in some way or the data type doesn't quite match
  //
  // the "new_owner" pointer has to have Py_DECREF called on it; it owns
  // the "new" array object created by PyArray_Cast
  //
  static PyArrayObject *tmp_arr;
  PyArrayObject *new_owner;
  tmp_arr = PyArray_GETCONTIGUOUS(array);
  new_owner = (PyArrayObject *) PyArray_Cast(tmp_arr, typenum);
  Py_DECREF(tmp_arr);
  return new_owner;
}


/* The PyInt variable is a PyLong in Python3.x.
 */
#if PY_MAJOR_VERSION >= 3
#define PyInt_AsLong PyLong_AsLong
#define PyInt_Check PyLong_Check
#endif

static PyObject *csolve(PyObject* self, PyObject *args, PyObject *kwargs)
{
  /* Expects a function call
   *     sol = csolve((m,n,p),c,Gx,Gi,Gp,h,dims,Ax,Ai,Ap,b)
   * where
   *
   * the triple (m,n,p) corresponds to:
   *    `m`: the rows of G
   *    `n`: the cols of G and A, must agree with the length of c
   *    `p`: the rows of A
   * `c` is a Numpy array of doubles
   * "G" is a sparse matrix in column compressed storage. "Gx" are the values,
   * "Gi" are the rows, and "Gp" are the column pointers.
   * `Gx` is a Numpy array of doubles
   * `Gi` is a Numpy array of ints
   * `Gp` is a Numpy array of ints
   * `h` is a Numpy array
   * `dims` is a dictionary with
   *    `dims['l']` an integer specifying the dimension of positive orthant cone
   *    `dims['q']` an *list* specifying dimensions of second-order cones
   *
   * "A" is an optional sparse matrix in column compressed storage. "Ax" are 
   * the values, "Ai" are the rows, and "Ap" are the column pointers.
   * `Ax` is a Numpy array of doubles
   * `Ai` is a Numpy array of ints
   * `Ap` is a Numpy array of ints
   * `b` is an optional argument, which is a Numpy array of doubles
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
  //matrix *c, *h, *b = NULL;
  //spmatrix *G, *A = NULL;
  
  PyArrayObject *Gx, *Gi, *Gp, *c, *h;
  PyArrayObject *Ax = NULL;
  PyArrayObject *Ai = NULL;
  PyArrayObject *Ap = NULL;
  PyArrayObject *b = NULL;
  PyObject *dims;
  idxint n;      // number or variables
  idxint m;      // number of conic variables
  idxint p = 0;  // number of equality constraints
  idxint ncones = 0; // number of cones
  idxint numConicVariables = 0;

  /* ECOS data structures */
  idxint l = 0;
  idxint *q = NULL;


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
  static char *kwlist[] = {"shape", "c", "Gx", "Gi", "Gp", "h", "dims", "Ax", "Ai", "Ap", "b", NULL};
  // parse the arguments and ensure they are the correct type
  // TODO: (ECHU) allow a "settings" struct
#ifdef DLONG
  static char *argparse_string = "(lll)O!O!O!O!O!O!|O!O!O!O!";
#else
  static char *argparse_string = "(iii)O!O!O!O!O!O!|O!O!O!O!";
#endif
    
  if( !PyArg_ParseTupleAndKeywords(args, kwargs, argparse_string, kwlist,
      &m, &n, &p,
      &PyArray_Type, &c,
      &PyArray_Type, &Gx,
      &PyArray_Type, &Gi,
      &PyArray_Type, &Gp,
      &PyArray_Type, &h,
      &PyDict_Type, &dims,
      &PyArray_Type, &Ax,
      &PyArray_Type, &Ai,
      &PyArray_Type, &Ap,
      &PyArray_Type, &b)
    ) { return NULL; }
  
  if (m < 0) {
    PyErr_SetString(PyExc_ValueError, "m must be a positive integer");
    return NULL;
  }

  if (n <= 0) {
    PyErr_SetString(PyExc_ValueError, "n must be a positive integer");
    return NULL;
  }
  
  if (p < 0) {
    PyErr_SetString(PyExc_ValueError, "p must be a positive integer");
    return NULL;
  }
  
  /* get the typenum for the primitive int and double types */
  int intType = getIntType();
  int doubleType = getDoubleType();

  /* set G */
  if( !PyArray_ISFLOAT(Gx) || PyArray_NDIM(Gx) != 1) {
    PyErr_SetString(PyExc_TypeError, "Gx must be a numpy array of floats");
    return NULL;
  }
  if( !PyArray_ISINTEGER(Gi) || PyArray_NDIM(Gi) != 1) {
    PyErr_SetString(PyExc_TypeError, "Gi must be a numpy array of ints");
    return NULL;
  }
  if( !PyArray_ISINTEGER(Gp) || PyArray_NDIM(Gp) != 1) {
    PyErr_SetString(PyExc_TypeError, "Gp must be a numpy array of ints");
    return NULL;
  }
  PyArrayObject *Gx_arr = getContiguous(Gx, doubleType);
  PyArrayObject *Gi_arr = getContiguous(Gi, intType);
  PyArrayObject *Gp_arr = getContiguous(Gp, intType);
  Gpr = (pfloat *) PyArray_DATA(Gx_arr);
  Gir = (idxint *) PyArray_DATA(Gi_arr);
  Gjc = (idxint *) PyArray_DATA(Gp_arr);

  /* set c */
  if (!PyArray_ISFLOAT(c) || PyArray_NDIM(c) != 1) {
      PyErr_SetString(PyExc_TypeError, "c must be a dense numpy array with one dimension");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      return NULL;
  }
  
  if (PyArray_DIM(c,0) != n){
      PyErr_SetString(PyExc_ValueError, "c has incompatible dimension with G");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      return NULL;
  }
  PyArrayObject *c_arr = getContiguous(c, doubleType);
  cpr = (pfloat *) PyArray_DATA(c_arr);

  /* set h */
  if (!PyArray_ISFLOAT(h) || PyArray_NDIM(h) != 1) {
      PyErr_SetString(PyExc_TypeError, "h must be a dense numpy array with one dimension");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr);
      return NULL;
  }


  if (PyArray_DIM(h,0) != m){
      PyErr_SetString(PyExc_ValueError, "h has incompatible dimension with G");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr);
      return NULL;
  }
  PyArrayObject *h_arr = getContiguous(h, doubleType);
  hpr = (pfloat *) PyArray_DATA(h_arr);

  /* get dims['l'] */
  PyObject *linearObj = PyDict_GetItemString(dims, "l");
  if(linearObj) {
    if(PyInt_Check(linearObj) && ((l = (idxint) PyInt_AsLong(linearObj)) >= 0)) {
        numConicVariables += l;
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['l'] ought to be a nonnegative integer");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
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
              numConicVariables += q[i];
          } else {
            PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list of positive integers");
            Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
            Py_DECREF(c_arr); Py_DECREF(h_arr);
            return NULL;
          }

      }
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
      return NULL;
    }
  }

  PyArrayObject *Ax_arr = NULL;
  PyArrayObject *Ai_arr = NULL;
  PyArrayObject *Ap_arr = NULL;
  PyArrayObject *b_arr = NULL;
  if(Ax && Ai && Ap && b) {
    /* set A */
    if( !PyArray_ISFLOAT(Ax) || PyArray_NDIM(Ax) != 1 ) {
      PyErr_SetString(PyExc_TypeError, "Ax must be a numpy array of floats");
      if(q) free(q);
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
      return NULL;
    }
    if( !PyArray_ISINTEGER(Ai) || PyArray_NDIM(Ai) != 1) {
      PyErr_SetString(PyExc_TypeError, "Ai must be a numpy array of ints");
      if(q) free(q);
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
      return NULL;
    }
    if( !PyArray_ISINTEGER(Ap) || PyArray_NDIM(Ap) != 1) {
      PyErr_SetString(PyExc_TypeError, "Ap must be a numpy array of ints");
      if(q) free(q);
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
      return NULL;
    }
    // if ((SpMatrix_Check(A) && SP_ID(A) != DOUBLE)){
    //     PyErr_SetString(PyExc_TypeError, "A must be a sparse 'd' matrix");
    //     if(q) free(q);
    //     Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
    //     Py_DECREF(c_arr); Py_DECREF(h_arr);
    //     return NULL;
    // }
    // if ((p = SP_NROWS(A)) < 0) {
    //     PyErr_SetString(PyExc_ValueError, "p must be a nonnegative integer");
    //     if(q) free(q);
    //     Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
    //     Py_DECREF(c_arr); Py_DECREF(h_arr);
    //     return NULL;
    // }
    // if (SP_NCOLS(A) != n) {
    //     PyErr_SetString(PyExc_ValueError, "A has incompatible dimension with c");
    //     if(q) free(q);
    //     Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
    //     Py_DECREF(c_arr); Py_DECREF(h_arr);
    //     return NULL;
    // }
    // if (p != 0) {
    //   Apr = SP_VALD(A);
    //   Air = SP_ROW(A);
    //   Ajc = SP_COL(A);
    // }
    Ax_arr = getContiguous(Ax, doubleType);
    Ai_arr = getContiguous(Ai, intType);
    Ap_arr = getContiguous(Ap, intType);
    Apr = (pfloat *) PyArray_DATA(Ax_arr);
    Air = (idxint *) PyArray_DATA(Ai_arr);
    Ajc = (idxint *) PyArray_DATA(Ap_arr);

    /* set b */
    // if (!Matrix_Check(b) || MAT_NCOLS(b) != 1 || MAT_ID(b) != DOUBLE) {
    //     PyErr_SetString(PyExc_TypeError, "b must be a dense 'd' matrix with one column");
    //     if(q) free(q);
    //     return NULL;
    // }
    // if (MAT_NROWS(b) != p){
    //     PyErr_SetString(PyExc_ValueError, "b has incompatible dimension with A");
    //     if(q) free(q);
    //     return NULL;
    // }
    // if (p != 0) {
    //   bpr = MAT_BUFD(b);
    // }
    if (!PyArray_ISFLOAT(b) || PyArray_NDIM(b) != 1) {
        PyErr_SetString(PyExc_TypeError, "b must be a dense numpy array with one dimension");
        if(q) free(q);
        Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
        Py_DECREF(c_arr); Py_DECREF(h_arr);
        Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
        return NULL;
    }
    if (PyArray_DIM(b,0) != p){
        PyErr_SetString(PyExc_ValueError, "b has incompatible dimension with A");
        if(q) free(q);
        Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
        Py_DECREF(c_arr); Py_DECREF(h_arr);
        Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
        return NULL;
    }
    b_arr = getContiguous(b, doubleType);
    bpr = (pfloat *) PyArray_DATA(b_arr);
  } else if (Ax || Ai || Ap || b) {
    // check that A and b are both supplied
    PyErr_SetString(PyExc_ValueError, "A and b arguments must be supplied together");
    if(q) free(q);
    Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
    Py_DECREF(c_arr); Py_DECREF(h_arr);
    return NULL;
  }


  /* check that sum(q) + l = m */
  if( numConicVariables != m ){
      PyErr_SetString(PyExc_ValueError, "Number of rows of G does not match dims.l+sum(dims.q)");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr); 
      if (b_arr) Py_DECREF(b_arr);
      if (Ax_arr) Py_DECREF(Ax_arr); 
      if (Ai_arr) Py_DECREF(Ai_arr); 
      if (Ap_arr) Py_DECREF(Ap_arr);
      return NULL;
  }
  
  /* This calls ECOS setup function. */
  mywork = ECOS_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, cpr, hpr, bpr);
  if( mywork == NULL ){
      PyErr_SetString(PyExc_RuntimeError, "Internal problem occurred in ECOS while setting up the problem.\nPlease send a bug report with data to Alexander Domahidi.\nEmail: domahidi@control.ee.ethz.ch");
      if(q) free(q);
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
      if (b_arr) Py_DECREF(b_arr);
      if (Ax_arr) Py_DECREF(Ax_arr); 
      if (Ai_arr) Py_DECREF(Ai_arr); 
      if (Ap_arr) Py_DECREF(Ap_arr);
      return NULL;
  }
  
  /* Solve! */
  idxint exitcode = ECOS_solve(mywork);

  /* create output (all data is *deep copied*) */
  // TODO: request CVXOPT API for constructing from existing pointer
  /* x */
  // matrix *x;
  // if(!(x = Matrix_New(n,1,DOUBLE)))
  //   return PyErr_NoMemory();
  // memcpy(MAT_BUFD(x), mywork->x, n*sizeof(double));
  npy_intp veclen[1];
  veclen[0] = n;
  PyObject *x = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->x);

  /* y */
  // matrix *y;
  // if(!(y = Matrix_New(p,1,DOUBLE)))
  //   return PyErr_NoMemory();
  // memcpy(MAT_BUFD(y), mywork->y, p*sizeof(double));
  veclen[0] = p;
  PyObject *y = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->y);
  
  
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
      case ECOS_NUMERICS:
          infostring = "Run into numerical problems";
          break;
      case ECOS_OUTCONE:
          infostring = "PROBLEM: Multipliers leaving the cone";
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

#if PROFILING > 0
  // give reference to infoDict
  Py_DECREF(tinfos);
#endif

  /* s */
  // matrix *s;
  // if(!(s = Matrix_New(m,1,DOUBLE)))
  //   return PyErr_NoMemory();
  // memcpy(MAT_BUFD(s), mywork->s, m*sizeof(double));
  veclen[0] = m;
  PyObject *s = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->s);
  
  /* z */
  // matrix *z;
  // if(!(z = Matrix_New(m,1,DOUBLE)))
  //   return PyErr_NoMemory();
  // memcpy(MAT_BUFD(z), mywork->z, m*sizeof(double));
  veclen[0] = m;
  PyObject *z = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->z);
  


  /* cleanup */
  ECOS_cleanup(mywork, 4);

  PyObject *returnDict = Py_BuildValue(
    "{s:O,s:O,s:O,s:O,s:O}",
    "x",x,
    "y",y,
    "z",z,
    "s",s,
    "info",infoDict);
  // give up ownership to the return dictionary
  Py_DECREF(x); Py_DECREF(y); Py_DECREF(z); Py_DECREF(s); Py_DECREF(infoDict);
  
  // no longer need pointers to arrays that held primitives
  if(q) free(q);
  Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
  Py_DECREF(c_arr); Py_DECREF(h_arr);
  if (b_arr) Py_DECREF(b_arr);
  if (Ax_arr) Py_DECREF(Ax_arr); 
  if (Ai_arr) Py_DECREF(Ai_arr); 
  if (Ap_arr) Py_DECREF(Ap_arr);

  return returnDict;
}

static PyMethodDef ECOSMethods[] =
{
  {"csolve", (PyCFunction)csolve, METH_VARARGS | METH_KEYWORDS,
    "Solve an SOCP using ECOS."},
  {NULL, NULL, 0, NULL} // sentinel
};

/* Module initialization */
#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef moduledef = {
          PyModuleDef_HEAD_INIT,
          "csolve",     /* m_name */
          "Solve an SOCP using ECOS.",  /* m_doc */
          -1,                  /* m_size */
          ECOSMethods,         /* m_methods */
          NULL,                /* m_reload */
          NULL,                /* m_traverse */
          NULL,                /* m_clear */
          NULL,                /* m_free */
  };
#endif

static PyObject* moduleinit(void)
{
  PyObject* m;

#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
#else
  m = Py_InitModule("_ecos", ECOSMethods);
#endif
  
  //if (import_array() < 0) return NULL; // for numpy arrays
  //if (import_cvxopt() < 0) return NULL; // for cvxopt support

  if (m == NULL)
    return NULL;

  return m;
};

#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_ecos(void)
  {
    import_array(); // for numpy arrays
    return moduleinit();
  }
#else
  PyMODINIT_FUNC init_ecos(void)
  {
    import_array(); // for numpy arrays
    moduleinit();
  }
#endif
