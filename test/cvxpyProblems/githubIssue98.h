#include "ecos.h"
#include "minunit.h"

static pfloat Gx[17] = {0.416757847405471, 2.136196095668454, 1.793435585194863, -1.,
  0.056266827226329, -1.640270808404989, 0.841747365656204, -1.,
  0.416757847405471, 2.136196095668454, 1.793435585194863, -1.,
  0.056266827226329, -1.640270808404989, 0.841747365656204, -1., -1.};
static idxint Gp[6] = {0, 4, 8, 12, 16, 17};
static idxint Gi[17] = {0, 1, 2, 7, 0, 1, 2, 8, 3, 4, 5, 9, 3, 4, 5, 10, 6};

static idxint q[1] = {5};
static pfloat c[5] = {0, 0, 0, 0, 1};
static pfloat h[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


static char * test_issue98()
{
    /**
     * minimize 0
     * s.t. 0 <= x <= 1
     */
    pwork *mywork;
    idxint exitflag;

    /* set up data */
    mywork = ECOS_setup(5, 11, 0,
        6, 1, q, 0,
        Gx, Gp, Gi,
        NULL, NULL, NULL,
        c, h, NULL);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork);
    }
    else exitflag = ECOS_FATAL;

    /* clean up memory */
    ECOS_cleanup(mywork, 0);

    mu_assert("githubIssue98-test: ECOS failed to produce outputflag OPTIMAL", exitflag == ECOS_OPTIMAL );
    return 0;
}
