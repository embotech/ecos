#include "ecos.h"
#include "minunit.h"
static char * test_infeasible2(){
idxint n = 3;
idxint m = 3;
idxint p = 2;
idxint l = 3;
idxint ncones = 0;
pfloat c[3] = {-1.0, 0.0, 0.0};
idxint Gjc[4] = {0, 0, 2, 4};
idxint Gir[4] = {0, 1, 0, 2};
pfloat Gpr[4] = {8.0, 1.0, 9.0, -1.0};
pfloat h[3] = {25.0, 0.0, 0.0};
idxint *q = NULL;
idxint Ajc[4] = {0,2,3,4};
idxint Air[4] = {0,1,1,1};
pfloat Apr[4] = {1.0, -1.0, 1.0, 1.0};
pfloat b[2] = {3.0, 0.0};

pwork *mywork;
idxint exitflag;

/* set up data */
#ifdef EXPCONE
mywork = ECOS_setup(n, m, p, l, ncones, q, 0, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
#else
mywork = ECOS_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
#endif

if( mywork != NULL ){
/* solve */
exitflag = ECOS_solve(mywork); }
else exitflag = ECOS_FATAL;

/* clean up memory */
ECOS_cleanup(mywork, 0);

mu_assert("infeasible2: ECOS failed to produce outputflag PRIMAL INFEASIBLE", exitflag == ECOS_PINF );
return 0;
}
