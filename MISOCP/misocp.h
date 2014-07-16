#ifndef __MISOCP_H__
#define __MISOCP_H__

#include "ecos.h"

// MISOCP configuration settings
#define INITIAL_NODES (100)
#define EXPAND_NODES (100)

// Preprocess flags
#define EMPTY (-2)
#define VALID (-1)

typedef struct misocp_node {
	// ret_code of -1 is a valid unsolved node
	// ret_code of -2 is an overwriteable node
	idxint ret_code;
	idxint last_split;
	
	pfloat lb;
	pfloat* constraints;
} misocp_node;

// Wrapper for mixed integer module
typedef struct misocp_pwork{
	// Mixed integer data
	idxint num_int;
	idxint num_nodes;
	misocp_node* nodes;
	misocp_node tmp_node;

	// ECOS data
	pwork* ecos_prob;
	
} misocp_pwork;


#endif