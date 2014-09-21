#include "ecos.h"

int main(){
	idxint n = 3;
	idxint m = 6;
	pfloat feas_Gx[8] = {1,-1,9999,-9999,1,1,1,1};
	idxint feas_Gp[4] = {0,4,6,8};
	idxint feas_Gi[8] = {0,1,2,3,2,4,3,5};

	pfloat feas_c[3] = {0,1,-1};
	pfloat feas_h[6] = {1,0,9999,0,10,12};

	// Answer: 
	pfloat x[3] = {0,10,0};

	idxint i, ret_code, pass;

	pwork* prob = ECOS_setup(
		n, m, 0, 
		m, 0, NULL,
		feas_Gx, feas_Gp, feas_Gi,
		NULL, NULL, NULL,
		feas_c, feas_h, NULL);

	ret_code = ECOS_solve(prob);
			
	ECOS_cleanup(prob,0);
	
	return 0;
}