#include "ecos.h"
#include "minunit.h"
static char * test_unboundedMaxSqrt(){
idxint n = 2;
idxint m = 3;
idxint p = 0;
idxint l = 0;
idxint ncones = 1;
pfloat c[2] = {0.0, -1.000000000000000000e+00};
idxint Gjc[3] = {0, 2, 3};
idxint Gir[3] = {0, 2, 1};
pfloat Gpr[3] = {-1.000000000000000000e+00, -1.000000000000000000e+00, -2.000000000000000000e+00};
pfloat h[3] = {1.000000000000000000e+00, 0.0, -1.000000000000000000e+00};
idxint q[1] = {3};
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
 
mu_assert("unboundedMaxSqrt: ECOS failed to produce outputflag DUAL INFEASIBLE", exitflag == ECOS_DINF );
return 0;
}
