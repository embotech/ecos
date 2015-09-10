#include "ecos.h"
#include "minunit.h"

/* Solves log(ax)-x for a = 0.3*/
idxint logax_x_n = 2;
idxint logax_x_m = 3;
idxint logax_x_p = 0;
idxint logax_x_l = 0;
idxint logax_x_ncones = 0;
idxint logax_x_nexc = 1;

pfloat logax_x_c[2] = {-1.0, 1.0};

idxint logax_x_Gjc[3] = {0,1,2};
idxint logax_x_Gir[2] = {0, 1};
pfloat logax_x_Gpr[2] = {-0.3, -1.0};

pfloat logax_x_h[3] = {0.0,0.0,1.0};

static char * test_log_ax_x(){
pwork *mywork;
idxint exitflag;
 
/* print test name */
printf("=================================== log(ax)-x ===================================\n");
 
/* set up data */
mywork = ECOS_setup(logax_x_n, logax_x_m, logax_x_p, logax_x_l, logax_x_ncones, NULL, logax_x_nexc,
                    logax_x_Gpr, logax_x_Gjc, logax_x_Gir,
                    NULL,NULL,NULL,
                    logax_x_c, logax_x_h, NULL);
if( mywork != NULL ){
/*set precision*/
	mywork->stgs->feastol = 1e-12;
/* solve */
exitflag = ECOS_solve(mywork); 
}
else exitflag = ECOS_FATAL;

/* Exact solution */
pfloat true_x = -log(0.3)/0.3;
pfloat ecos_x = mywork->x[0];
pfloat abs_err = fabs(true_x-ecos_x);

/* clean up memory */
ECOS_cleanup(mywork, 0);

mu_assert("log_ax_x: ECOS failed to produce outputflag OPTIMAL", exitflag == ECOS_OPTIMAL );
mu_assert("log_ax_x: ECOS failed to produce the desired precision", abs_err < 1e-11);

return 0;
}
