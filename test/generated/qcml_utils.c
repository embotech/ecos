#include <stdlib.h>
#include "qcml_utils.h"

/* free a qc_socp structure */
qc_socp *qc_socp_free(qc_socp *data)
{
  if(data) {
    if (data->q) free(data->q);
    if (data->Gx) free(data->Gx);
    if (data->Gp) free(data->Gp);
    if (data->Gi) free(data->Gi);
    if (data->Ax) free(data->Ax);
    if (data->Ap) free(data->Ap);
    if (data->Ai) free(data->Ai);
    if (data->c) free(data->c);
    if (data->h) free(data->h);
    if (data->b) free(data->b);
    free(data);
  }
  return NULL;
}


/*
 * private utility functions
 * duplicates (some) SuiteSparse csparse functionality
 */

/* free a coo matrix */
qc_matrix *qc_spfree(qc_matrix *A)
{
  if(A) {
    if (A->i) free(A->i);
    if (A->j) free(A->j);
    if (A->v) free(A->v);
    free(A);
  }
  return NULL;
}

/* allocate memory for a coo matrix */
qc_matrix *qc_spalloc (idxint m, idxint n, idxint nnz, int triplet)
{
  /* allocate the coo struct */
  qc_matrix *A = (qc_matrix *) calloc (1, sizeof (qc_matrix)) ;
  if (!A) return (NULL) ;                 /* out of memory */
  A->m = m ;                              /* define dimensions and nzmax */
  A->n = n ;
  A->nnz = triplet ? nnz : -1 ;
  A->j = (idxint *) malloc (triplet ? (nnz * sizeof(idxint)) : ((n+1) * sizeof(idxint))) ;
  A->i = (idxint *) malloc (nnz * sizeof (idxint)) ;
  A->v = (pfloat *) malloc (nnz * sizeof (pfloat)) ;
  return ((!A->v || !A->i || !A->j) ? qc_spfree(A) : A);
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
void cumsum (idxint *p, idxint *c, idxint n)
{
  /* performs a cumulative sum; may overflow if exceed idxint storage (4GB
   * worth of nonzerors)
   *
   * although could write c[i] += c[i-1], we don't. the extra workspace allows
   * the compiler to optimize the code
   */
  idxint i, nz = 0 ;
  if (!p || !c) return ;          /* check inputs */
  for (i = 0 ; i < n ; i++)
  {
      p [i] = nz ;
      nz += c [i] ;
      c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
  }
  p [n] = nz ;
}

/* remove duplicate entries from A. doesn't realloc memory.
 * should only be called on coo matrix that is moonlighting as CSC matrix.
 * not publicly exposed, so it doesn't really matter. 
 */
int remove_dup (qc_matrix *A)
{
  idxint i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w ;
  pfloat *Ax ;
  m = A->m ; n = A->n ; Ap = A->j ; Ai = A->i ; Ax = A->v ;
  w = (idxint *) malloc (m * sizeof (idxint)) ;   /* get workspace */
  if (!w) return 0;                           /* out of memory */
  for (i = 0 ; i < m ; i++) w [i] = -1 ;      /* row i not yet seen */
  for (j = 0 ; j < n ; j++)
  {
    q = nz ;                                /* column j will start at q */
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      i = Ai [p] ;                        /* A(i,j) is nonzero */
      if (w [i] >= q)
      {
        Ax [w [i]] += Ax [p] ;          /* A(i,j) is a duplicate */
      }
      else
      {
        w [i] = nz ;                    /* record where row i occurs */
        Ai [nz] = i ;                   /* keep A(i,j) */
        Ax [nz++] = Ax [p] ;
      }
    }
    Ap [j] = q ;                            /* record start of column j */
  }
  Ap [n] = nz ;                               /* finalize A */
  A->nnz = nz ;
  free (w) ;                                  /* free workspace */
  return 1 ;                                  /* success! */
}

/* compresses i,j,v pointers into CSC format but stores it into a coo struct
 * removes duplicate entires as a last step */
qc_matrix *qc_compress (const qc_matrix *T)
{
  idxint m, n, nnz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
  pfloat *Cx, *Tx ;
  qc_matrix *C ;
  
  if (T->v == NULL) return NULL;
  
  m = T->m ; n = T->n ; Ti = T->i ; Tj = T->j ; Tx = T->v ; nnz = T->nnz ;
  C = qc_spalloc(m, n, nnz, QC_CSC) ;    /* allocate memory for CSC format */
  if (!C) return NULL;
  
  /* create temporary workspace */
  w = (idxint *) calloc (n, sizeof(idxint));
  if (!w) {
    free(C->i);
    free(C->j);
    free(C->v);
    free(C);
    return NULL;  
  }                         
  Cp = C->j ; Ci = C->i ; Cx = C->v ;
  for (k = 0 ; k < nnz ; k++) w [Tj [k]]++ ;          /* column counts */
  cumsum (Cp, w, n) ;                                 /* column pointers */

  for (k = 0 ; k < nnz ; k++)
  {
      Ci [p = w [Tj [k]]++] = Ti [k] ;   /* A(i,j) is the pth entry in C */
      if (Cx) Cx [p] = Tx [k] ;
  }
  
  /* success; free w and return C */
  free(w);
  if (remove_dup(C)) return C;    
  return NULL;  
}
