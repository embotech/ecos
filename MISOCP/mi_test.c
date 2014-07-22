#include "misocp.h"
#include "misocp.c"
#include "misocp_preproc.c"
#include "ecos.h"
#include "timer.h"

static idxint n = 2;
static idxint m = 4;
static pfloat feas_Gx[6] = {2.0, 3.0, -1.0, 1.0, 1.0, 4.0};
static idxint feas_Gp[3] = {0, 4, 6};
static idxint feas_Gi[6] = {0, 1, 2, 3, 0, 1};

static pfloat feas_c[2] = {-1, -1};
static pfloat feas_h[4] = {4, 12, 0 , 1};

int test_one(){
	idxint i, ret_code;
	
	misocp_pwork* prob = misocp_setup(
		n, m, 0, 
    	m, 0, NULL,
    	feas_Gx, feas_Gp, feas_Gi,
    	NULL, NULL, NULL,
    	feas_c, feas_h, NULL, 2);

	timer tfactor;
	tic(&tfactor);
	ret_code = misocp_solve(prob);
	pfloat t = toc(&tfactor);

	printf ("\nFound soln: ");
	for (i=0; i<n; ++i) printf ("%f," , prob->best_x[i]); printf ("\n");

	misocp_cleanup(prob);
			
	printf("Return code was: %u, Terminated in: %f seconds\n", ret_code,t);
	return ret_code;
}


int main(){
	printf("Test 1: %u\n", test_one());

	/*pwork* prob = ECOS_setup(
		n, m, 0, 
    	m, 0, NULL,
    	feas_Gx, feas_Gp, feas_Gi,
    	NULL, NULL, NULL,
    	feas_c, feas_h, NULL);

	//unset_equilibration(prob);
	prob->h[2] = -1.0;
	prob->h[3] = 1.0;
	//set_equilibration(prob);	
	ECOS_solve(prob);
	ECOS_solve(prob);
	ECOS_solve(prob);
	ECOS_solve(prob);
	ECOS_solve(prob);
	printf ("\nFound soln: ");
	for (i=0; i<n; ++i) printf ("%f," , prob->x[i]); printf ("\n");*/
	
	return 0;
}
