#include "misocp.h"
#include "misocp.c"
#include "misocp_preproc.c"
#include "ecos.h"
#include "timer.h"

int test_1(){
	idxint n = 2;
	idxint m = 4;
	pfloat feas_Gx[6] = {2.0, 3.0, -1.0, 1.0, 1.0, 4.0};
	idxint feas_Gp[3] = {0, 4, 6};
	idxint feas_Gi[6] = {0, 1, 2, 3, 0, 1};

	pfloat feas_c[2] = {-1, -1};
	pfloat feas_h[4] = {4, 12, 0 , 1};

	idxint i, ret_code, pass;
	
	misocp_pwork* prob = misocp_setup(
		n, m, 0, 
    	m, 0, NULL,
    	feas_Gx, feas_Gp, feas_Gi,
    	NULL, NULL, NULL,
    	feas_c, feas_h, NULL, 1);

	ret_code = misocp_solve(prob);
	
	pass = 1;
	// Answer is x = [1, 2]
	pfloat x[2] = {1.0, 2.0};
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
	}
	
	return pass;
}

int test_2(){
	idxint n = 2;
	idxint m = 4;
	pfloat feas_Gx[6] = {2.0, 3.0, -1.0, 1.0, 1.0, 4.0};
	idxint feas_Gp[3] = {0, 4, 6};
	idxint feas_Gi[6] = {0, 1, 2, 3, 0, 1};

	pfloat feas_c[2] = {-1., -1.};
	pfloat feas_h[4] = {4., 12., 0. , 1.};

	idxint i, ret_code, pass;

	misocp_pwork* prob = misocp_setup(
		n, m, 0, 
    	m, 0, NULL,
    	feas_Gx, feas_Gp, feas_Gi,
    	NULL, NULL, NULL,
    	feas_c, feas_h, NULL, 2);

	ret_code = misocp_solve(prob);
	
	pass = 1;
	// Answer is x = [1, 1]
	pfloat x[2] = {1.0, 1.0};
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
	}

	return pass;
}


int main(){
	printf("Pass: 1, Fail: 0\n");
	printf("Test 1: %u\n", test_1());
	printf("Test 2: %u\n", test_2());	
	
	return 0;
}
