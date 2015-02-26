#include "ecos.h"
#include "minunit.h"
#include "norm/norm.h"
#include "quad_over_lin/quad_over_lin.h"
#include "sq_norm/sq_norm.h"
#include "sum_sq/sum_sq.h"
#include "inv_pos/inv_pos.h"
#include "qcml_utils.h"

static char * test_norm()
{

    pwork *mywork;
    idxint exitflag;
    qc_socp *data;

    data = qc_norm2socp(NULL, NULL);

    /* set up data */
    mywork = ECOS_setup(data->n, data->m, data->p,
        data->l, data->nsoc, data->q, 0,
        data->Gx, data->Gp, data->Gi,
        data->Ax, data->Ap, data->Ai,
        data->c, data->h, data->b);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork);
    }
    else exitflag = ECOS_FATAL;

    /* clean up memory */
    ECOS_cleanup(mywork, 0);
    qc_socp_free(data);

    mu_assert("quadratics-norm-test: ECOS failed to produce output flag OPTIMAL", exitflag == ECOS_OPTIMAL );
    return 0;
}

static char * test_quad_over_lin()
{

    pwork *mywork;
    idxint exitflag;
    qc_socp *data;

    data = qc_quad_over_lin2socp(NULL, NULL);

    /* set up data */
    mywork = ECOS_setup(data->n, data->m, data->p,
        data->l, data->nsoc, data->q, 0,
        data->Gx, data->Gp, data->Gi,
        data->Ax, data->Ap, data->Ai,
        data->c, data->h, data->b);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork);
    }
    else exitflag = ECOS_FATAL;

    /* clean up memory */
    ECOS_cleanup(mywork, 0);
    qc_socp_free(data);

    mu_assert("quadratics-quad-over-lin-test: ECOS failed to produce output flag OPTIMAL", exitflag == ECOS_OPTIMAL );
    return 0;
}

static char * test_sq_norm()
{

    pwork *mywork;
    idxint exitflag;
    qc_socp *data;

    data = qc_sq_norm2socp(NULL, NULL);

    /* set up data */
    mywork = ECOS_setup(data->n, data->m, data->p,
        data->l, data->nsoc, data->q, 0,
        data->Gx, data->Gp, data->Gi,
        data->Ax, data->Ap, data->Ai,
        data->c, data->h, data->b);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork);
    }
    else exitflag = ECOS_FATAL;

    /* clean up memory */
    ECOS_cleanup(mywork, 0);
    qc_socp_free(data);

    mu_assert("quadratics-sq-norm-test: ECOS failed to produce output flag OPTIMAL", exitflag == ECOS_OPTIMAL );
    return 0;
}

static char * test_sum_sq()
{

    pwork *mywork;
    idxint exitflag;
    qc_socp *data;

    data = qc_sum_sq2socp(NULL, NULL);

    /* set up data */
    mywork = ECOS_setup(data->n, data->m, data->p,
        data->l, data->nsoc, data->q, 0,
        data->Gx, data->Gp, data->Gi,
        data->Ax, data->Ap, data->Ai,
        data->c, data->h, data->b);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork);
    }
    else exitflag = ECOS_FATAL;

    /* clean up memory */
    ECOS_cleanup(mywork, 0);
    qc_socp_free(data);

    mu_assert("quadratics-sum-sq-test: ECOS failed to produce output flag OPTIMAL", exitflag == ECOS_OPTIMAL );
    return 0;
}

static char * test_inv_pos()
{

    pwork *mywork;
    idxint exitflag;
    qc_socp *data;

    data = qc_inv_pos2socp(NULL, NULL);

    /* set up data */
    mywork = ECOS_setup(data->n, data->m, data->p,
        data->l, data->nsoc, data->q, 0, 
        data->Gx, data->Gp, data->Gi,
        data->Ax, data->Ap, data->Ai,
        data->c, data->h, data->b);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork);
    }
    else exitflag = ECOS_FATAL;

    /* clean up memory */
    ECOS_cleanup(mywork, 0);
    qc_socp_free(data);

    mu_assert("inv-pos-test: ECOS failed to produce at least something close to OPTIMAL", exitflag == ECOS_OPTIMAL+ECOS_INACC_OFFSET || exitflag == ECOS_OPTIMAL );
    return 0;
}

