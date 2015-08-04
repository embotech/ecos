#include "ecos.h"
#include "minunit.h"
static char * test_unboundedLP1(){
idxint n = 2;
idxint m = 4;
idxint p = 0;
idxint l = 4;
idxint ncones = 0;
pfloat c[2] = {-1.000000000000000000e+00, -1.000000000000000000e+00};
idxint Gjc[3] = {0, 3, 6};
idxint Gir[6] = {0, 1, 2, 0, 1, 3};
pfloat Gpr[6] = {-2.000000000000000000e+00, -1.000000000000000222e+00, -1.000000000000000000e+00, -9.999999999999998890e-01, -3.000000000000000000e+00, -1.000000000000000000e+00};
pfloat h[4] = {-1.000000000000000000e+00, -1.000000000000000000e+00, 0.0, 0.0};
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
 
mu_assert("unboundedLP1: ECOS failed to produce outputflag UNBOUNDED", exitflag == ECOS_DINF );
return 0;
}
