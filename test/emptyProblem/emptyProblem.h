#include "ecos.h"
#include "minunit.h"
static char * test_emptyProblem(){
idxint n = 0;
idxint m = 0;
idxint p = 0;
idxint l = 0;
idxint ncones = 0;
pfloat *c = NULL;
idxint *Gjc = NULL;
idxint *Gir = NULL;
pfloat *Gpr = NULL;
pfloat *h = NULL;
idxint *q = NULL;
idxint *Ajc = NULL;
idxint *Air = NULL;
pfloat *Apr = NULL;
pfloat *b = NULL;
 
pwork *mywork;
idxint exitflag;
 
/* set up data */
mywork = ECOS_setup(n, m, p, l, ncones, q, 0,Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
if( mywork != NULL ){
/* solve */
exitflag = ECOS_solve(mywork); }
else exitflag = ECOS_FATAL;
 
/* clean up memory */
ECOS_cleanup(mywork, 0);
 
mu_assert("emptyProblem: ECOS failed to produce outputflag OPTIMAL", exitflag == ECOS_OPTIMAL );
return 0;
}
