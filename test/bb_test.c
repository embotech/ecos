#include "timer.h"
#include "ecos.h"
#include "ecos_bb.h"

int test_basic(){
	idxint n = 2;
	idxint m = 1;
	pfloat feas_Ax[1] = {1.0};
	idxint feas_Ap[3] = {0, 1, 1};
	idxint feas_Ai[1] = {0};

	pfloat feas_Gx[1] = {-1.0};
	idxint feas_Gp[3] = {0, 0, 1};
	idxint feas_Gi[1] = {0};

	pfloat feas_c[2] = {1., 1.};
	pfloat feas_b[1] = {0.};
	pfloat feas_h[1] = {-0.5};

	idxint bool_idx[1] = {1};

	/* Answer: */
	pfloat x[2] = {0.0, 0.0};

	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ECOS_BB_setup(
                                      n, m, 1,
                                      1, 0, NULL, 0,
                                      feas_Gx, feas_Gp, feas_Gi,
                                      feas_Ax, feas_Ap, feas_Ai,
                                      feas_c, feas_h, feas_b, 1, bool_idx, 0 , NULL, NULL);

	printf("Passed setup \n");

	prob->stgs->verbose = 0;
	ret_code = ECOS_BB_solve(prob);

	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i], prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	ECOS_BB_cleanup(prob, 0);
	return pass;
}

int test_1a_bool(){
	idxint n = 2;
	idxint m = 2;
	pfloat feas_Gx[4] = {2.0, 3.0, 1.0, 4.0};
	idxint feas_Gp[3] = {0, 2, 4};
	idxint feas_Gi[4] = {0, 1, 0, 1};

	pfloat feas_c[2] = {-1.1, -1.};
	pfloat feas_h[2] = {4., 12.};

	idxint bool_idx[1] = {0};

	/* Answer: */
	pfloat x[2] = {1.0, 2.0};

	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0,
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 1, bool_idx, 0 , NULL, NULL);

	printf("Passed setup \n");

	prob->stgs->verbose = 0;
	ret_code = ECOS_BB_solve(prob);

	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i], prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	ECOS_BB_cleanup(prob, 0);
	return pass;
}

int test_1a_int(){
	idxint n = 2;
	idxint m = 4;
	pfloat feas_Gx[6] = {2.0, 3.0, -1, 1.0, 4.0, -1};
	idxint feas_Gp[3] = {0, 3, 6};
	idxint feas_Gi[6] = {0, 1, 2, 0, 1, 3};

	pfloat feas_c[2] = {-1., -1.1};
	pfloat feas_h[4] = {4., 12., 0. , 0.};

	idxint int_idx[2] = {0,1};

	/* Answer: */
	pfloat x[2] = {0.0, 3.0};

	idxint i, ret_code, pass;
	
	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 0 , NULL, 2, int_idx, NULL);
	prob->stgs->verbose = 0;
	printf("Passed setup \n");

	ret_code = ECOS_BB_solve(prob);
	
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	ECOS_BB_cleanup(prob, 0);
	return pass;
}

int test_1b(){
	idxint n = 2;
	idxint m = 2;
	pfloat feas_Gx[4] = {2.0, 3.0, 1.0, 4.0};
	idxint feas_Gp[3] = {0, 2, 4};
	idxint feas_Gi[4] = {0, 1, 0, 1};

	pfloat feas_c[2] = {-1., -1.};
	pfloat feas_h[2] = {4., 12.};

	idxint bool_idx[1] = {1};

	/* Answer: */
	pfloat x[2] = {1.5, 1.0};

	idxint i, ret_code, pass;
	
	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 1, bool_idx, 0 , NULL, NULL);
	prob->stgs->verbose = 0;
	printf("Passed setup \n");

	ret_code = ECOS_BB_solve(prob);
	
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	ECOS_BB_cleanup(prob, 0);
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

	idxint bool_idx[2] = {0,1};

	/* Answer: */
	pfloat x[2] = {1.0, 1.0};
	
	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 2, bool_idx, 0 , NULL, NULL);
	prob->stgs->verbose = 0;
	ret_code = ECOS_BB_solve(prob);
	
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	ECOS_BB_cleanup(prob, 0);
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

	idxint bool_idx[1] = {0};

	/* Answer: */
	pfloat x[2] = {0.0, 3.0};

	idxint i, ret_code, pass;
	
	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 1, bool_idx, 0 , NULL, NULL);
	prob->stgs->verbose = 0;
	ret_code = ECOS_BB_solve(prob);
	
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	ECOS_BB_cleanup(prob, 0);
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

	idxint bool_idx[6] = {0,1,2,3,4,5};

	/* Answer: */ 
	pfloat x[6] = {0,1,1,0,0,0};

	idxint i, ret_code, pass;

	
	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 6, bool_idx, 0, NULL, NULL);
	prob->stgs->verbose = 0;
	ret_code = ECOS_BB_solve(prob);

	pass = 1;
	
	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	
	ECOS_BB_cleanup(prob, 0);
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

	idxint bool_idx[6] = {0,1,2,3,4,5};

	/* Answer: */
	pfloat x[6] = {0,0,1,1,1,0};

	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
	feas_c, feas_h, NULL, 6, bool_idx, 0, NULL, NULL);
	prob->stgs->verbose = 0;
	ret_code = ECOS_BB_solve(prob);
	
	pass = 1;
	
	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}

	ECOS_BB_cleanup(prob, 0);

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

	idxint bool_idx[6] = {0,1,2,3,4,5};

	/* Answer: */
	pfloat x[6] = {0,1,1,1,1,0};

	idxint i, ret_code, pass;

	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 6, bool_idx, 0, NULL, NULL);
	prob->stgs->verbose = 0;
	printf("Obj:");
	for (i=0; i<n; ++i){
		printf("%f ", prob->ecos_prob->c[i]);
	}
	printf("\n");

	ret_code = ECOS_BB_solve(prob);
	
	printf("Obj:");
	for (i=0; i<n; ++i){
		printf("%f ", prob->ecos_prob->c[i]);
	}
	printf("\n");
		
	pass = 1;

	printf("Soln:");
	for (i=0; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	ECOS_BB_cleanup(prob, 0);
	return pass;
}

int test_7(){
	idxint n = 15;
	idxint m = 29;
	pfloat feas_Gx[120] = {9999,-9999,9999,-9999,9999,-9999,9999,-9999,9999,-9999,-3.5008,3.5008,-0.4504,0.4504,-0.8764999999999999,0.8764999999999999,-0.1088,0.1088,1,1,-1,-8.4095,8.4095,-1.0107,1.0107,-1.686,1.686,-0.3525,0.3525,1,1,-1,-15.1987,15.1987,-2.0203,2.0203,-2.3932,2.3932,-0.6233,0.6233,1,1,-1,-22.5405,22.5405,-3.1862,3.1862,-2.8749,2.8749,-0.7923,0.7923,1,1,-1,-29.2639,29.2639,-4.3096,4.3096,-3.0189,3.0189,-0.8116,0.8116,1,1,-1,3.5008,-3.5008,0.4504,-0.4504,0.8764999999999999,-0.8764999999999999,0.1088,-0.1088,1,1,-1,8.4095,-8.4095,1.0107,-1.0107,1.686,-1.686,0.3525,-0.3525,1,1,-1,15.1987,-15.1987,2.0203,-2.0203,2.3932,-2.3932,0.6233,-0.6233,1,1,-1,22.5405,-22.5405,3.1862,-3.1862,2.8749,-2.8749,0.7923,-0.7923,1,1,-1,29.2639,-29.2639,4.3096,-4.3096,3.0189,-3.0189,0.8116,-0.8116,1,1,-1};
	idxint feas_Gp[16] = {0,2,4,6,8,10,21,32,43,54,65,76,87,98,109,120};
	idxint feas_Gi[120] = {8,9,10,11,12,13,14,15,16,17,0,1,2,3,4,5,6,7,8,18,19,0,1,2,3,4,5,6,7,10,18,20,0,1,2,3,4,5,6,7,12,18,21,0,1,2,3,4,5,6,7,14,18,22,0,1,2,3,4,5,6,7,16,18,23,0,1,2,3,4,5,6,7,9,18,24,0,1,2,3,4,5,6,7,11,18,25,0,1,2,3,4,5,6,7,13,18,26,0,1,2,3,4,5,6,7,15,18,27,0,1,2,3,4,5,6,7,17,18,28};

	pfloat feas_c[15] = {0,0,0,0,0,0.127,0.9134,0.6324,0.0975,0.2785,0.873,0.0866,0.3676,0.9025,0.7215};
	pfloat feas_h[29] = {-729.9349999999999,789.9349999999999,-71.015,131.015,-89.66,149.66,-1.165,61.165,9999,0,9999,0,9999,0,9999,0,9999,0,150,0,0,0,0,0,0,0,0,0,0};

	idxint bool_idx[5] = {0,1,2,3,4};

	/* Answer: */ 
	pfloat x[15] = {0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,32.383266,0.00,0.00,0.00,
		0.00,0.00,0.00};
	pfloat x2[15] = {0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,86.798858,
		0.000000,0.000000,0.000000};

	idxint i, ret_code, pass;

	timer t;

	ecos_bb_pwork* prob = ECOS_BB_setup(
		n, m, 0, 
		m, 0, NULL, 0,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL, 5, bool_idx, 0, NULL, NULL);
	prob->stgs->verbose = 0;
	tic(&t);
	ret_code = ECOS_BB_solve(prob);
	pfloat msRuntime = toc(&t);

	pass = 1;

	printf("Soln: ");
	for (i=5; i<n; ++i){
		pass &= float_eqls(x[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	printf("\nRuntime: %f\n", msRuntime);
	
	updateDataEntry_h(prob, 0, 789.935);
	updateDataEntry_h(prob, 1, -729.935);
	updateDataEntry_h(prob, 2, 131.015);
	updateDataEntry_h(prob, 3, -71.015);
	updateDataEntry_h(prob, 4, 149.66);
	updateDataEntry_h(prob, 5, -89.66);
	updateDataEntry_h(prob, 6, 61.165);
	updateDataEntry_h(prob, 7, -1.165);

	tic(&t);
	ret_code = ECOS_BB_solve(prob);
	msRuntime = toc(&t);

	printf("Soln2: ");
	for (i=5; i<n; ++i){
		pass &= float_eqls(x2[i] ,prob->x[i], prob->stgs->integer_tol );
		printf("%f ", prob->x[i]);
	}
	printf("\nRuntime: %f\n", msRuntime);
	
	ECOS_BB_cleanup(prob, 0);

	return pass;
}


int main(){

	printf("Pass: 1, Fail: 0\n=============\n");
	printf("\n\nTest 1aBoolean: %u\n=============\n", test_1a_bool());
	printf("\n\nTest 1aInteger: %u\n=============\n", test_1a_int());
	printf("\n\nTest 1b: %u\n=============\n", test_1b());
	printf("\n\nTest 2: %u\n=============\n", test_2());	
	printf("\n\nTest 3: %u\n=============\n", test_3());	
	printf("\n\nTest 4: %u\n=============\n", test_4());
	printf("\n\nTest 5: %u\n=============\n", test_5());
	printf("\n\nTest 6: %u\n=============\n", test_6());	
	printf("\n\nTest 7: %u\n=============\n", test_7());			
	printf("\n\nTest basic: %u\n=============\n", test_basic());			
	return 0;
}
