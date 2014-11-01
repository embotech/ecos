/* Check that we are clean against numpy 1.7 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include "ecos.h"
#include "ecos_bb.h"
#include "numpy/arrayobject.h"
/*
 * Define INLINE for MSVC compatibility.
 */
#ifdef _MSC_VER
  #define INLINE __inline
#else
  #define INLINE inline
#endif

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
/* #include "cvxopt.h" */

/* ECHU: Note, Python3.x may require special handling for the int and double
 * types. */
static INLINE int getIntType(void) {
  switch(sizeof(idxint)) {
    case 1: return NPY_INT8;
    case 2: return NPY_INT16;
    case 4: return NPY_INT32;
    case 8: return NPY_INT64;
    default: return NPY_INT32;  /* defaults to 4 byte int */
  }
}

static INLINE int getDoubleType(void) {
  /* ECHU: known bug, if pfloat isn't "double", will cause aliasing in memory */
  return NPY_DOUBLE;
}

static INLINE PyArrayObject *getContiguous(PyArrayObject *array, int typenum) {
  /* gets the pointer to the block of contiguous C memory
   * the overhead should be small unless the numpy array has been
   * reordered in some way or the data type doesn't quite match
   *
   * the "new_owner" pointer has to have Py_DECREF called on it; it owns
   * the "new" array object created by PyArray_Cast
   */
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

static PyObject *version(PyObject* self)
{
  return Py_BuildValue("s",ECOS_VERSION);
}

static int checkNonnegativeInt(const char *key, idxint val) {
    if (val >= 0) {
      return 0;
    }
    PyErr_Format(PyExc_ValueError, "'%s' must be a nonnegative integer", key);
    return -1;
}

static int checkPositiveFloat(const char *key, pfloat val) {
    if (val > 0) {
      return 0;
    }
    PyErr_Format(PyExc_ValueError, "'%s' must be a positive float", key);
    return -1;
}

static PyObject *csolve(PyObject* self, PyObject *args, PyObject *kwargs)
{
  /* Expects a function call
   *     sol = csolve((m,n,p),c,Gx,Gi,Gp,h,dims,Ax,Ai,Ap,b,**kwargs)
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
   * other optional arguments are:
   *     `feastol`: the tolerance on the primal and dual residual
   *     `abstol`: the absolute tolerance on the duality gap
   *     `reltol`: the relative tolerance on the duality gap
   *     `feastol_inacc`: the tolerance on the primal and dual residual if reduced precisions
   *     `abstol_inacc`: the absolute tolerance on the duality gap if reduced precision
   *     `reltolL_inacc`: the relative tolerance on the duality gap if reduced precision
   *     `max_iters`: the maximum numer of iterations.
   *     `verbose`: signals to print on non zero value.
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
  /* ECHU: below is for CVXOPT
   * matrix *c, *h, *b = NULL;
   * spmatrix *G, *A = NULL;
   */

  /* BEGIN VARIABLE DECLARATIONS */
  PyArrayObject *Gx, *Gi, *Gp, *c, *h;
  PyArrayObject *bool_idx = NULL;
  PyArrayObject *Ax = NULL;
  PyArrayObject *Ai = NULL;
  PyArrayObject *Ap = NULL;
  PyArrayObject *b = NULL;
  PyObject *dims = NULL;
  PyObject *verbose = NULL;
  idxint n;           /* number or variables            */
  idxint m;           /* number of conic variables      */
  idxint p = 0;       /* number of equality constraints */
  idxint ncones = 0;  /* number of cones                */
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

  idxint *bool_vars_idx = NULL;
  idxint num_bool = 0;

  /* Default ECOS settings */
  settings opts_ecos;

  pwork* mywork;
  ecos_bb_pwork* myecos_bb_work = NULL;

  idxint i;
  static char *kwlist[] = {"shape", "c", "Gx", "Gi", "Gp", "h", "dims",
      "Ax", "Ai", "Ap", "b",
      "verbose", "feastol", "abstol", "reltol",
      "feastol_inacc", "abstol_inacc", "reltol_inacc",
      "max_iters", "integer_vars_idx", "num_bool", NULL};
  int intType, doubleType;

  /* parse the arguments and ensure they are the correct type */
#ifdef DLONG
  static char *argparse_string = "(lll)O!O!O!O!O!O!|O!O!O!O!O!ddddddlO!i";
#else
  static char *argparse_string = "(iii)O!O!O!O!O!O!|O!O!O!O!O!ddddddiO!i";
#endif
  PyArrayObject *Gx_arr, *Gi_arr, *Gp_arr;
  PyArrayObject *c_arr;
  PyArrayObject *h_arr;
  PyObject *linearObj;
  PyObject *socObj;
  PyArrayObject *Ax_arr = NULL;
  PyArrayObject *Ai_arr = NULL;
  PyArrayObject *Ap_arr = NULL;
  PyArrayObject *b_arr = NULL;
  PyArrayObject *bool_idx_arr = NULL;

  idxint exitcode, numerr = 0;
  npy_intp veclen[1];
  PyObject *x, *y, *z, *s;
  const char* infostring;
  PyObject *infoDict, *tinfos;
  PyObject *returnDict = NULL;
  /* END VARIABLE DECLARATIONS */

  /* Default ECOS settings */
  opts_ecos.feastol = FEASTOL;
  opts_ecos.reltol = RELTOL;
  opts_ecos.abstol = ABSTOL;
  opts_ecos.feastol_inacc = FTOL_INACC;
  opts_ecos.abstol_inacc = ATOL_INACC;
  opts_ecos.reltol_inacc = RTOL_INACC;
  opts_ecos.maxit = MAXIT;
  opts_ecos.verbose = VERBOSE;

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
      &PyArray_Type, &b,
      &PyBool_Type, &verbose,
      &opts_ecos.feastol,
      &opts_ecos.abstol,
      &opts_ecos.reltol,
      &opts_ecos.feastol_inacc,
      &opts_ecos.abstol_inacc,
      &opts_ecos.reltol_inacc,
      &opts_ecos.maxit,
      &PyArray_Type, &bool_idx)
    ) { return NULL; }

  if (checkNonnegativeInt("m", m) < 0) return NULL;
  if (checkNonnegativeInt("n", n) < 0) return NULL;
  if (checkNonnegativeInt("p", p) < 0) return NULL;

  /* check the opts*/
  if (verbose)
      opts_ecos.verbose = (idxint) PyObject_IsTrue(verbose);
  if (checkNonnegativeInt("maxit", opts_ecos.maxit) < 0) return NULL;
  if (checkPositiveFloat("abstol", opts_ecos.abstol) < 0) return NULL;
  if (checkPositiveFloat("feastol", opts_ecos.feastol) < 0) return NULL;
  if (checkPositiveFloat("reltol", opts_ecos.reltol) < 0) return NULL;
  if (checkPositiveFloat("abstol_inacc", opts_ecos.abstol_inacc) < 0) return NULL;
  if (checkPositiveFloat("feastol_inacc", opts_ecos.feastol_inacc) < 0) return NULL;
  if (checkPositiveFloat("reltol_inacc", opts_ecos.reltol_inacc) < 0) return NULL;

  /* get the typenum for the primitive int and double types */
  intType = getIntType();
  doubleType = getDoubleType();

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
  Gx_arr = getContiguous(Gx, doubleType);
  Gi_arr = getContiguous(Gi, intType);
  Gp_arr = getContiguous(Gp, intType);
  Gpr = (pfloat *) PyArray_DATA(Gx_arr);
  Gir = (idxint *) PyArray_DATA(Gi_arr);
  Gjc = (idxint *) PyArray_DATA(Gp_arr);

  /* set c */
  if (!PyArray_ISFLOAT(c) || PyArray_NDIM(c) != 1) {
      PyErr_SetString(PyExc_TypeError, "c must be a dense numpy float array with one dimension");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      return NULL;
  }

  if (PyArray_DIM(c,0) != n){
      PyErr_SetString(PyExc_ValueError, "c has incompatible dimension with G");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      return NULL;
  }
  c_arr = getContiguous(c, doubleType);
  cpr = (pfloat *) PyArray_DATA(c_arr);

  /* set h */
  if (!PyArray_ISFLOAT(h) || PyArray_NDIM(h) != 1) {
      PyErr_SetString(PyExc_TypeError, "h must be a dense numpy float array with one dimension");
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
  h_arr = getContiguous(h, doubleType);
  hpr = (pfloat *) PyArray_DATA(h_arr);

  /* get dims['l'] */
  linearObj = PyDict_GetItemString(dims, "l");
  if(linearObj) {
    if ( (PyInt_Check(linearObj) && ((l = (idxint) PyInt_AsLong(linearObj)) >= 0)) ||
         (PyLong_Check(linearObj) && ((l = PyLong_AsLong(linearObj)) >= 0)) ){
        numConicVariables += l;
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['l'] ought to be a nonnegative integer");
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
      return NULL;
    }
  }

  /* get dims['q'] */
  socObj = PyDict_GetItemString(dims, "q");
  if(socObj) {
    if (PyList_Check(socObj)) {
      ncones = PyList_Size(socObj);
      q = calloc(ncones, sizeof(idxint));
      for (i = 0; i < ncones; ++i) {
          PyObject *qi = PyList_GetItem(socObj, i);
          if( (PyInt_Check(qi) && ((q[i] = (idxint) PyInt_AsLong(qi)) > 0)) ||
              (PyLong_Check(qi) && ((q[i] = PyLong_AsLong(qi)) > 0)) ) {
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
    /* if ((SpMatrix_Check(A) && SP_ID(A) != DOUBLE)){
     *     PyErr_SetString(PyExc_TypeError, "A must be a sparse 'd' matrix");
     *     if(q) free(q);
     *     Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
     *     Py_DECREF(c_arr); Py_DECREF(h_arr);
     *     return NULL;
     * }
     * if ((p = SP_NROWS(A)) < 0) {
     *     PyErr_SetString(PyExc_ValueError, "p must be a nonnegative integer");
     *     if(q) free(q);
     *     Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
     *     Py_DECREF(c_arr); Py_DECREF(h_arr);
     *     return NULL;
     * }
     * if (SP_NCOLS(A) != n) {
     *     PyErr_SetString(PyExc_ValueError, "A has incompatible dimension with c");
     *     if(q) free(q);
     *     Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
     *     Py_DECREF(c_arr); Py_DECREF(h_arr);
     *     return NULL;
     * }
     * if (p != 0) {
     *   Apr = SP_VALD(A);
     *   Air = SP_ROW(A);
     *   Ajc = SP_COL(A);
     * }
     */
    Ax_arr = getContiguous(Ax, doubleType);
    Ai_arr = getContiguous(Ai, intType);
    Ap_arr = getContiguous(Ap, intType);
    Apr = (pfloat *) PyArray_DATA(Ax_arr);
    Air = (idxint *) PyArray_DATA(Ai_arr);
    Ajc = (idxint *) PyArray_DATA(Ap_arr);

    /* set b */
    /* if (!Matrix_Check(b) || MAT_NCOLS(b) != 1 || MAT_ID(b) != DOUBLE) {
     *     PyErr_SetString(PyExc_TypeError, "b must be a dense 'd' matrix with one column");
     *     if(q) free(q);
     *     return NULL;
     * }
     * if (MAT_NROWS(b) != p){
     *     PyErr_SetString(PyExc_ValueError, "b has incompatible dimension with A");
     *     if(q) free(q);
     *     return NULL;
     * }
     * if (p != 0) {
     *   bpr = MAT_BUFD(b);
     * }
     */
    if (!PyArray_ISFLOAT(b) || PyArray_NDIM(b) != 1) {
        PyErr_SetString(PyExc_TypeError, "b must be a dense numpy float array with one dimension");
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
    /* check that A and b are both supplied */
    PyErr_SetString(PyExc_ValueError, "A and b arguments must be supplied together");
    if(q) free(q);
    Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
    Py_DECREF(c_arr); Py_DECREF(h_arr);
    return NULL;
  }

  /* check that sum(q) + l = m */
  if( numConicVariables != m ){
      PyErr_SetString(PyExc_ValueError, "Number of rows of G does not match dims.l+sum(dims.q)");
      if (q) free(q);
      Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
      Py_DECREF(c_arr); Py_DECREF(h_arr);
      if (b_arr) Py_DECREF(b_arr);
      if (Ax_arr) Py_DECREF(Ax_arr);
      if (Ai_arr) Py_DECREF(Ai_arr);
      if (Ap_arr) Py_DECREF(Ap_arr);
      return NULL;
  }

  if (bool_idx){
    num_bool = (idxint) PyArray_DIM(bool_idx,0);
    bool_idx_arr = getContiguous(bool_idx, intType);
    bool_vars_idx = (idxint *) PyArray_DATA(bool_idx_arr);

    /* This calls ECOS setup function. */
    myecos_bb_work = ecos_bb_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, cpr, hpr, bpr, num_bool, bool_vars_idx);
    if( myecos_bb_work == NULL ){
        PyErr_SetString(PyExc_RuntimeError, "Internal problem occurred in ECOS_BB while setting up the problem.\nPlease send a bug report with data to Alexander Domahidi.\nEmail: domahidi@control.ee.ethz.ch");
        if(q) free(q);
        Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
        Py_DECREF(c_arr); Py_DECREF(h_arr);
        if (b_arr) Py_DECREF(b_arr);
        if (Ax_arr) Py_DECREF(Ax_arr);
        if (Ai_arr) Py_DECREF(Ai_arr);
        if (Ap_arr) Py_DECREF(Ap_arr);
        Py_DECREF(bool_idx);
        return NULL;
    }

    mywork = myecos_bb_work->ecos_prob;

    /* Set settings for ECOS. */
    mywork->stgs->verbose = opts_ecos.verbose;
    mywork->stgs->abstol = opts_ecos.abstol;
    mywork->stgs->feastol = opts_ecos.feastol;
    mywork->stgs->reltol = opts_ecos.reltol;
    mywork->stgs->abstol_inacc = opts_ecos.abstol_inacc;
    mywork->stgs->feastol_inacc = opts_ecos.feastol_inacc;
    mywork->stgs->reltol_inacc = opts_ecos.reltol_inacc;
    mywork->stgs->maxit = opts_ecos.maxit;

    /* Solve! */
    exitcode = ecos_bb_solve(myecos_bb_work);

  } else{
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

    /* Set settings for ECOS. */
    mywork->stgs->verbose = opts_ecos.verbose;
    mywork->stgs->abstol = opts_ecos.abstol;
    mywork->stgs->feastol = opts_ecos.feastol;
    mywork->stgs->reltol = opts_ecos.reltol;
    mywork->stgs->abstol_inacc = opts_ecos.abstol_inacc;
    mywork->stgs->feastol_inacc = opts_ecos.feastol_inacc;
    mywork->stgs->reltol_inacc = opts_ecos.reltol_inacc;
    mywork->stgs->maxit = opts_ecos.maxit;

    /* Solve! */
    exitcode = ECOS_solve(mywork);
  }

  /* create output (all data is *deep copied*) */
  /* TODO: request CVXOPT API for constructing from existing pointer */
  /* x */
  /* matrix *x;
   * if(!(x = Matrix_New(n,1,DOUBLE)))
   *   return PyErr_NoMemory();
   * memcpy(MAT_BUFD(x), mywork->x, n*sizeof(double));
   */
  veclen[0] = n;
  x = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->x);

  /* y */
  /* matrix *y;
   * if(!(y = Matrix_New(p,1,DOUBLE)))
   *   return PyErr_NoMemory();
   * memcpy(MAT_BUFD(y), mywork->y, p*sizeof(double));
   */
  veclen[0] = p;
  y = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->y);

  if (bool_idx){
    /* info dict */
    /* infostring */
    switch( exitcode ){
        case MI_OPTIMAL_SOLN:
            infostring = "Optimal branch and bound solution found";
            break;
        case MI_MAXITER_FEASIBLE_SOLN:
            infostring = "Maximum iterations reached with feasible solution found";
            break;
        case MI_MAXITER_NO_SOLN:
            infostring = "Maximum iterations reached with no feasible solution found";
            break;
        case MI_INFEASIBLE:
            infostring = "Problem is infeasible";
            break;
        default:
            infostring = "UNKNOWN PROBLEM IN BRANCH AND BOUND SOLVER";
    }
  } else {
    /* info dict */
    /* infostring */
    switch( exitcode ){
        case ECOS_OPTIMAL:
            infostring = "Optimal solution found";
            break;
        case ECOS_OPTIMAL + ECOS_INACC_OFFSET:
            infostring = "Close to optimal solution found";
            break;
        case ECOS_MAXIT:
            infostring = "Maximum number of iterations reached";
            break;
        case ECOS_PINF:
            infostring = "Primal infeasible";
            break;
        case ECOS_PINF + ECOS_INACC_OFFSET:
            infostring = "Close to primal infeasible";
            break;
        case ECOS_DINF:
            infostring = "Dual infeasible";
            break;
        case ECOS_DINF + ECOS_INACC_OFFSET:
            infostring = "Close to dual infeasible";
            break;
        case ECOS_NUMERICS:
            infostring = "Run into numerical problems";
            break;
        case ECOS_OUTCONE:
            infostring = "PROBLEM: Multipliers leaving the cone";
            break;
        case ECOS_FATAL:
            infostring = "PROBLEM: Fatal error during initialization";
            break;
        default:
            infostring = "UNKNOWN PROBLEM IN SOLVER";
    }

    /* numerical errors */
    if( (exitcode == ECOS_NUMERICS) || (exitcode == ECOS_OUTCONE) || (exitcode == ECOS_FATAL) ){
        numerr = 1;
    } 
  }

  /* timings */
#if PROFILING > 0
	tinfos = Py_BuildValue(
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

  infoDict = Py_BuildValue(
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
  /* give reference to infoDict */
  Py_DECREF(tinfos);
#endif

  /* s */
  /* matrix *s;
   * if(!(s = Matrix_New(m,1,DOUBLE)))
   *   return PyErr_NoMemory();
   * memcpy(MAT_BUFD(s), mywork->s, m*sizeof(double));
   */
  veclen[0] = m;
  s = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->s);

  /* z */
  /* matrix *z;
   * if(!(z = Matrix_New(m,1,DOUBLE)))
   *   return PyErr_NoMemory();
   * memcpy(MAT_BUFD(z), mywork->z, m*sizeof(double));
   */
  veclen[0] = m;
  z = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, mywork->z);

  /* cleanup */
  if (bool_idx){
    ecos_bb_cleanup(myecos_bb_work, 4);
  }else{
    ECOS_cleanup(mywork, 4);  
  }
  
  returnDict = Py_BuildValue(
    "{s:O,s:O,s:O,s:O,s:O}",
    "x",x,
    "y",y,
    "z",z,
    "s",s,
    "info",infoDict);
  /* give up ownership to the return dictionary */
  Py_DECREF(x); Py_DECREF(y); Py_DECREF(z); Py_DECREF(s); Py_DECREF(infoDict);

  /* no longer need pointers to arrays that held primitives */
  if(q) free(q);
  Py_DECREF(Gx_arr); Py_DECREF(Gi_arr); Py_DECREF(Gp_arr);
  Py_DECREF(c_arr); Py_DECREF(h_arr);
  if (b_arr) Py_DECREF(b_arr);
  if (Ax_arr) Py_DECREF(Ax_arr);
  if (Ai_arr) Py_DECREF(Ai_arr);
  if (Ap_arr) Py_DECREF(Ap_arr);
  if (bool_idx_arr) Py_DECREF(bool_idx_arr);

  return returnDict;
}

static PyMethodDef ECOSMethods[] =
{
  {"csolve", (PyCFunction)csolve, METH_VARARGS | METH_KEYWORDS,
    "Solve an SOCP using ECOS."},
  {"version", (PyCFunction)version, METH_NOARGS, "Version number for ECOS."},
  {NULL, NULL, 0, NULL} /* sentinel */
};

/* Module initialization */
#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef moduledef = {
          PyModuleDef_HEAD_INIT,
          "_ecos",             /* m_name */
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

  /*if (import_array() < 0) return NULL;  */ /* for numpy arrays */
  /*if (import_cvxopt() < 0) return NULL; */ /* for cvxopt support */

  if (m == NULL)
    return NULL;

  return m;
};

#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit__ecos(void)
  {
    import_array(); /* for numpy arrays */
    return moduleinit();
  }
#else
  PyMODINIT_FUNC init_ecos(void)
  {
    import_array(); /* for numpy arrays */
    moduleinit();
  }
#endif
