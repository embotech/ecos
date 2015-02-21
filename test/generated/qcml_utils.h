/*
 * TODO: MIT/BSD license this code
 *
 * This is an ANSI C compatible file that contains the basic utility and data
 * structures needed for QCML.
*/
#ifndef __QCML_UTILS_H__
#define __QCML_UTILS_H__

#include "glblopts.h"

#define QC_CSC 0
#define QC_COO 1

#ifdef __cplusplus
extern "C" {
#endif

/* external coordinate (i,j,v) format for sparse matrices
 * TODO: if i, j are NULL but v is not, the matrix is assumed to be 
 * dense in C (row-major) ordering
 */
typedef struct coo {
  pfloat *v;  /* nonzero values */
  idxint *i;    /* row pointer    */
  idxint *j;    /* col pointer    */
  idxint nnz;   /* number of nonzeros, -1 if CSC */
  idxint m;     /* number of rows in the matrix  */
  idxint n;     /* number of cols in the matrix  */
} qc_matrix;

/* the socp data struct */
typedef struct socp {
  idxint n;     /* number of variables             */
  idxint m;     /* number of cone constraints      */
  idxint p;     /* number of equality constraints  */
  idxint l;     /* number of linear cones          */
  idxint nsoc;  /* number of second-order cones    */
  idxint *q;    /* list of second-order cone sizes */
  pfloat *Gx; /* nonzero values of G (in CSC)    */
  idxint *Gp;   /* column pointers of G (in CSC)   */
  idxint *Gi;   /* row values of G (in CSC)        */
  pfloat *Ax; /* nonzero values of A (in CSC)    */
  idxint *Ap;   /* column pointers of A (in CSC)   */
  idxint *Ai;   /* row values of A (in CSC)        */
  pfloat *c;  /* c vector (dense)                */
  pfloat *h;  /* h vector (dense)                */
  pfloat *b;  /* b vector (dense)                */
} qc_socp;

/* free an allocated socp data struct
 *     socp is allocated by corresponding prob2socp function
 */
qc_socp *qc_socp_free(qc_socp *data);

/* free an allocated sparse matrix */
qc_matrix *qc_spfree(qc_matrix *A);

/* compress a COO matrix into CSC, allocates new matrix memory */
qc_matrix *qc_compress(const qc_matrix *T);

#ifdef __cplusplus
}
#endif
  
#endif /* __QCML_UTILS_H__ */