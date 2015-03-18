#include "ecos.h"
#include "minunit.h"

static pfloat feas_Gx[2] = {1, -1};
static idxint feas_Gp[2] = {0, 2};
static idxint feas_Gi[2] = {0, 1};

static pfloat feas_c[1] = {0};
static pfloat feas_h[2] = {1, 0};

static char * test_feas()
{
    /**
     * minimize 0
     * s.t. 0 <= x <= 1
     */
    pwork *mywork;
    idxint exitflag;

    /* set up data */
    mywork = ECOS_setup(1, 2, 0,
        2, 0, NULL, 0,
        feas_Gx, feas_Gp, feas_Gi,
        NULL, NULL, NULL,
        feas_c, feas_h, NULL);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork); 
    }
    else exitflag = ECOS_FATAL;
        
    /* clean up memory */
    ECOS_cleanup(mywork, 0);
    
    mu_assert("feas-test: ECOS failed to produce outputflag OPTIMAL", exitflag == ECOS_OPTIMAL );
    return 0;
}
