#include "ecos.h"
#include "minunit.h"
static char * test_infeasible1(){
idxint n = 1;
idxint m = 2;
idxint p = 0;
idxint l = 2;
idxint ncones = 0;
pfloat c[1] = {-1.000000000000000000e+00};
idxint Gjc[2] = {0, 2};
idxint Gir[2] = {0, 1};
pfloat Gpr[2] = {-1.000000000000000000e+00, 1.000000000000000000e+00};
pfloat h[2] = {-2.000000000000000000e+00, 1.000000000000000000e+00};
idxint *q = NULL;
idxint *Ajc = NULL;
idxint *Air = NULL;
pfloat *Apr = NULL;
pfloat *b = NULL;
 
pwork *mywork;
idxint exitflag;
 
/* set up data */
mywork = ECOS_setup(n, m, p, l, ncones, q, 0, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
if( mywork != NULL ){
/* solve */
exitflag = ECOS_solve(mywork); }
else exitflag = ECOS_FATAL;
 
/* clean up memory */
ECOS_cleanup(mywork, 0);
 
mu_assert("infeasible1: ECOS failed to produce outputflag PRIMAL INFEASIBLE", exitflag == ECOS_PINF );
return 0;
}
