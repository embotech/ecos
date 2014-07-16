#ifndef __MISOCP_H__
#define __MISOCP_H__

typedef struct misocp_node {
	idxint ret_code;
	idxint can_overwrite;
	pfloat up_bound;
	pfloat * constraints;
} misocp_node;

// Wrapper for mixed integer module
typedef struct misocp_pwork{
	// Mixed integer data
	idxint num_int;
	idxint num_nodes;
	misocp_node * nodes;

	// ECOS data
	pwork * ecos_prob;
	
} misocp_pwork;


#endif