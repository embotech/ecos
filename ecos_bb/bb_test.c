#include "timer.h"
#include "ecos_bb.h"
#include "ecos_bb.c"
#include "ecos_bb_preproc.c"

int test_1(){
	idxint n = 2;
	idxint m = 2;
	pfloat feas_Gx[4] = {2.0, 3.0, 1.0, 4.0};
	idxint feas_Gp[3] = {0, 2, 4};
	idxint feas_Gi[4] = {0, 1, 0, 1};

	pfloat feas_c[2] = {-1., -1.};
	pfloat feas_h[2] = {4., 12.};

	// Answer:
	pfloat x[2] = {1.0, 2.0};

	idxint i, ret_code, pass;
	
	ecos_bb_pwork* prob = ecos_bb_setup(
		n, m, 0, 
		m, 0, NULL,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 1);

	ret_code = ecos_bb_solve(prob);
	
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
		printf("%f,", prob->best_x[i]);
	}
	ecos_bb_cleanup(prob);
	return pass;
}

int test_2(){
	idxint n = 2;
	idxint m = 2;
	pfloat feas_Gx[4] = {2.0, 3.0, 1.0, 4.0};
	idxint feas_Gp[3] = {0, 2, 4};
	idxint feas_Gi[4] = {0, 1, 0, 1};

	pfloat feas_c[2] = {-1., -1.};
	pfloat feas_h[2] = {4., 12.};

	// Answer:
	pfloat x[2] = {1.0, 1.0};
	
	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ecos_bb_setup(
		n, m, 0, 
		m, 0, NULL,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 2);

	ret_code = ecos_bb_solve(prob);
	
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
		printf("%f,", prob->best_x[i]);
	}
	ecos_bb_cleanup(prob);
	return pass;
}


int test_3(){
	idxint n = 2;
	idxint m = 2;
	pfloat feas_Gx[4] = {2.0, 3.0, 1.0, 4.0};
	idxint feas_Gp[3] = {0, 2, 4};
	idxint feas_Gi[4] = {0, 1, 0, 1};

	pfloat feas_c[2] = {1., -1.};
	pfloat feas_h[2] = {4., 12.};

	// Answer:
	pfloat x[2] = {0.0, 3.0};

	idxint i, ret_code, pass;
	
	ecos_bb_pwork* prob = ecos_bb_setup(
		n, m, 0, 
		m, 0, NULL,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 1);

	ret_code = ecos_bb_solve(prob);
	
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
		printf("%f,", prob->best_x[i]);
	}
	ecos_bb_cleanup(prob);
	return pass;
}

int test_4(){
	idxint n = 6;
	idxint m = 3;
	pfloat feas_Gx[18] = {2,5,-5,-6,3,1,3,-1,-4,-4,-3,2,-1,2,-2,2,-1,1};
	idxint feas_Gp[7] = {0,3,6,9,12,15,18};
	idxint feas_Gi[18] = {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};

	pfloat feas_c[6] = {3, 5, 6, 9, 10, 10};
	pfloat feas_h[3] = {-2, 2, -3};

	// Answer: 
	pfloat x[6] = {0,1,1,0,0,0};

	idxint i, ret_code, pass;

	
	ecos_bb_pwork* prob = ecos_bb_setup(
		n, m, 0, 
		m, 0, NULL,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 6);

	ret_code = ecos_bb_solve(prob);

	pass = 1;
	
	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
		printf("%f,", prob->best_x[i]);
	}
	
	ecos_bb_cleanup(prob);
	return pass;
}


int test_5(){
	idxint n = 6;
	idxint m = 3;
	pfloat feas_Gx[18] = {2,5,-5,-6,3,1,3,-1,-4,-4,-3,2,-1,2,-2,2,-1,1};
	idxint feas_Gp[7] = {0,3,6,9,12,15,18};
	idxint feas_Gi[18] = {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};

	pfloat feas_c[6] = {3, 5, 6, 9, -10, -10};
	pfloat feas_h[3] = {-2, 2, -3};

	// Answer: 
	pfloat x[6] = {0,0,1,1,1,0};

	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ecos_bb_setup(
		n, m, 0, 
		m, 0, NULL,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
	feas_c, feas_h, NULL, 6);

	ret_code = ecos_bb_solve(prob);
	
	pass = 1;
	
	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
		printf("%f,", prob->best_x[i]);
	}

	ecos_bb_cleanup(prob);

	return pass;
}

int test_6(){
	idxint n = 6;
	idxint m = 3;
	pfloat feas_Gx[18] = {2,5,-5,-6,3,1,3,-1,-4,-4,-3,2,-1,2,-2,2,-1,1};
	idxint feas_Gp[7] = {0,3,6,9,12,15,18};
	idxint feas_Gi[18] = {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};

	pfloat feas_c[6] = {-3, -5, 6, 9, 10, -10};
	pfloat feas_h[3] = {-2, 1, -3};

	// Answer: 
	pfloat x[6] = {0,1,1,1,1,0};

	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ecos_bb_setup(
		n, m, 0, 
		m, 0, NULL,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 6);

	ret_code = ecos_bb_solve(prob);
	
	pass = 1;
	
	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->best_x[i]);
		printf("%f,", prob->best_x[i]);
	}
	ecos_bb_cleanup(prob);
	return pass;
}

int main(){

	printf("Pass: 1, Fail: 0\n=============\n");
	printf("\nTest 1: %u\n=============\n", test_1());
	printf("\nTest 2: %u\n=============\n", test_2());	
	printf("\nTest 3: %u\n=============\n", test_3());	
	printf("\nTest 4: %u\n=============\n", test_4());
	printf("\nTest 5: %u\n=============\n", test_5());
	printf("\nTest 6: %u\n=============\n", test_6());		
	
	return 0;
}
