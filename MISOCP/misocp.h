#ifndef __MISOCP_H__
#define __MISOCP_H__

#include "ecos.h"

// MISOCP configuration settings
#define MI_ABS_EPS (1E-7)
#define MI_REL_EPS (1E-2)
#define MI_MAXITER (1000)
#define MI_INT_TOL (1E-5)

// Flags
#define MI_SOLVED_NON_BRANCHABLE (3)
#define MI_SOLVED_BRANCHABLE (2)
#define MI_NOT_SOLVED (1)
#define MI_FREE (0)

#define MI_ONE (1)
#define MI_ZERO (0)
#define MI_STAR (-1)

#define INFINITY (1E999)

typedef struct node {
	char status;
	pfloat L;
	pfloat U;
	idxint split_idx;
} misocp_node;

// Wrapper for mixed integer module
typedef struct misocp_pwork{
	// Mixed integer data
	idxint num_int_vars;
	
	node* nodes;
	char* node_ids;
	
	pfloat* best_x;
    pfloat global_U;
    pfloat global_L;

	// ECOS data
	pwork* ecos_prob;

	// Tmp data
	char* tmp_node_id;
	idxint iter;
	
} misocp_pwork;

// Calculate the offset into the node_id array
inline pfloat* get_node_id(idxint idx, misocp_pwork* prob){
    return &prob->node_ids[prob->num_int_vars * idx];
}

inline pfloat round(pfloat number){
    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

inline pfloat abs(pfloat a, pfloat b){
	return a < b ? (b-a) : (a-b);
}

inline idxint float_eqls(pfloat a, pfloat b){
	return abs(a - b) < MI_INT_TOL;
}

inline pfloat min(pfloat a, pfloat b){
	return a < b ? a : b;
}



#endif