#include <stdlib.h>
#include "sum_sq.h"

/* ----------------------- BEGIN GENERATED CODE --------------------------- */
qc_socp * qc_sum_sq2socp()
{
    /*
     * maps 'params' into the C socp data type
     * 'params' ought to contain:
     */

    /* all local variables */
    idxint i;  /* loop index */
    idxint *q_ptr;
    idxint *A_row_ptr, *A_col_ptr;
    idxint *G_row_ptr, *G_col_ptr;
    pfloat *A_data_ptr, *G_data_ptr;
    idxint nnzA, nnzG;
    qc_matrix *G_csc;
    qc_matrix G_coo; 

    /* allocate socp data structure */
    qc_socp * data = (qc_socp *) calloc(1, sizeof(qc_socp));
    if (!data) return qc_socp_free(data);

    /* allocate problem dimensions */
    data->p = 0;
    data->m = 12;
    data->n = 6;

    /* allocate c vector */
    data->c = (pfloat *) calloc(data->n, sizeof(pfloat));
    if (!data->c) return qc_socp_free(data);

    /* allocate h vector */
    data->h = (pfloat *) calloc(data->m, sizeof(pfloat));
    if (!data->h) return qc_socp_free(data);

    /* allocate b vector */
    data->b = NULL;

    /* allocate G matrix */
    nnzG =  + 12;
    data->Gx = (pfloat *) malloc(nnzG * sizeof(pfloat));
    data->Gp = (idxint *) malloc(nnzG * sizeof(idxint));
    data->Gi = (idxint *) malloc(nnzG * sizeof(idxint));
    if ((!data->Gx) || (!data->Gp) || (!data->Gi)) return qc_socp_free(data);
    G_data_ptr = data->Gx;
    G_row_ptr = data->Gi;
    G_col_ptr = data->Gp;

    /* allocate A matrix */
    nnzA = 0;
    data->Ax = NULL;
    data->Ap = NULL;
    data->Ai = NULL;
    A_data_ptr = data->Ax;
    A_row_ptr = data->Ai;
    A_col_ptr = data->Ap;

    /* allocate the cone sizes */
    data->l = 3;
    data->nsoc = 3;
    data->q = (idxint *) malloc(data->nsoc * sizeof(idxint));
    if(!data->q) return qc_socp_free(data);

    /* initialize the cone */
    q_ptr = data->q;
    for(i = 0; i < 3; ++i) *q_ptr++ = 3;

    /* stuffing the objective vector */
    for(i = 0; i < 3; ++i) data->c[i + 0] = 1;

    /* for the constraint -1*x <= 0 */
    for(i = 0; i < 3; ++i) *G_row_ptr++ = i;
    for(i = 3; i < 6; ++i) *G_col_ptr++ = i;
    for(i = 0; i < 3; ++i) *G_data_ptr++ = -1;

    /* for the SOC product constraint norm(1 + -1*_t0, 2.0*x) <= 1 + _t0 */
    for(i = 5; i < 14; i+=3) *G_row_ptr++ = i;
    for(i = 3; i < 6; ++i) *G_col_ptr++ = i;
    for(i = 0; i < 3; ++i) *G_data_ptr++ = -2.0;
    for(i = 0; i < 3; ++i) data->h[3 * i + 4] = 1;
    for(i = 4; i < 13; i+=3) *G_row_ptr++ = i;
    for(i = 0; i < 3; ++i) *G_col_ptr++ = i;
    for(i = 0; i < 3; ++i) *G_data_ptr++ = 1;
    for(i = 0; i < 3; ++i) data->h[3 * i + 3] = 1;
    for(i = 3; i < 12; i+=3) *G_row_ptr++ = i;
    for(i = 0; i < 3; ++i) *G_col_ptr++ = i;
    for(i = 0; i < 3; ++i) *G_data_ptr++ = -1;

    /* convert G and A ptrs into a qc_matrix */
    G_coo.m = data->m; G_coo.n = data->n; G_coo.nnz = nnzG;
    G_coo.i = data->Gi;
    G_coo.j = data->Gp;
    G_coo.v = data->Gx;

    /* convert the matrices to column compressed form */
    G_csc = qc_compress(&G_coo);
    if (!G_csc) return qc_socp_free(data);
    /* free memory used for COO matrix, so it can be reassigned later */
    free(data->Gi);
    free(data->Gp);
    free(data->Gx);
    /* reassign into data, pointer now owned by data */
    data->Gi = G_csc->i;
    data->Gp = G_csc->j;
    data->Gx = G_csc->v;
    /* only free temp CSC pointer, but not its data */
    free(G_csc);

    return data;
}

void qc_socp2sum_sq(pfloat * x, sum_sq_vars * vars, const sum_sq_dims * dims)
{
    /*
     * recovers the problem variables from the solver variable 'x'
     * assumes the variables struct is externally allocated
     * the user must keep track of the variable length;
     */
    vars->x = x + 3;  /* length 3 */
}
/* ------------------------ END GENERATED CODE ---------------------------- */
