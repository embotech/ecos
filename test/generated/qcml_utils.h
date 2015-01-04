/*
 * TODO: MIT/BSD license this code
 *
 * This is an ANSI C compatible file that contains the basic utility and data
 * structures needed for QCML.
*/
#ifndef __QCML_UTILS_H__
#define __QCML_UTILS_H__

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
  double *v;  /* nonzero values */
  long *i;    /* row pointer    */
  long *j;    /* col pointer    */
  long nnz;   /* number of nonzeros, -1 if CSC */
  long m;     /* number of rows in the matrix  */
  long n;     /* number of cols in the matrix  */
} qc_matrix;

/* the socp data struct */
typedef struct socp {
  long n;     /* number of variables             */
  long m;     /* number of cone constraints      */
  long p;     /* number of equality constraints  */
  long l;     /* number of linear cones          */
  long nsoc;  /* number of second-order cones    */
  long *q;    /* list of second-order cone sizes */
  double *Gx; /* nonzero values of G (in CSC)    */
  long *Gp;   /* column pointers of G (in CSC)   */
  long *Gi;   /* row values of G (in CSC)        */
  double *Ax; /* nonzero values of A (in CSC)    */
  long *Ap;   /* column pointers of A (in CSC)   */
  long *Ai;   /* row values of A (in CSC)        */
  double *c;  /* c vector (dense)                */
  double *h;  /* h vector (dense)                */
  double *b;  /* b vector (dense)                */
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