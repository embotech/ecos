#ifndef __ecos_bb_H__
#define __ecos_bb_H__

#include "ecos.h"
#include "spla.h"
#include "glblopts.h"

/* Print verbosity */
#define MI_PRINTLEVEL (1)

/* ecos_bb configuration settings */
#define MI_ABS_EPS (1E-7)
#define MI_REL_EPS (1E-2)
#define MI_MAXITER (100)
#define MI_INT_TOL (1E-5)

/* Flags */
#define MI_SOLVED_NON_BRANCHABLE (3)
#define MI_SOLVED_BRANCHABLE (2)
#define MI_NOT_SOLVED (1)
#define MI_FREE (0)

#define MI_ONE (1)
#define MI_ZERO (0)
#define MI_STAR (-1)

#define MI_OPTIMAL_SOLN (0)
#define MI_MAXITER_FEASIBLE_SOLN (-1)
#define MI_MAXITER_NO_SOLN (1)
#define MI_INFEASIBLE (2)

#define MAX_FLOAT_INT (8388608) //Max integer and all smaller integer representable by single precision

#ifdef __cplusplus
extern "C" {
#endif

//#define INFINITY (1.0/0.0)

typedef struct node {
	char status;
	pfloat L;
	pfloat U;
	idxint split_idx;
	pfloat split_val;
} node;

/* Wrapper for mixed integer module */
typedef struct ecos_bb_pwork{
	/* Mixed integer data */
	idxint num_bool_vars;
	idxint num_int_vars;
	idxint maxiter;
	
	node* nodes;
	char* bool_node_ids;
	pfloat* int_node_ids;

	idxint* bool_vars_idx;
	idxint* int_vars_idx;

	/* ECOS data */
	pwork* ecos_prob;

	/* Modified pointers to ecos internals */
	/* Use these to edit or reset the h variables */
	spmat* A;
	spmat* G;
	pfloat* c;
	pfloat* b;
	pfloat* h;

	/* best iterate seen so far */
    /* variables */
    pfloat* x;  /* primal variables                    */
    pfloat* y;  /* multipliers for equality constaints */
    pfloat* z;  /* multipliers for conic inequalities  */
    pfloat* s;  /* slacks for conic inequalities       */
    pfloat kap; /* kappa (homogeneous embedding)       */
	pfloat tau; /* tau (homogeneous embedding)         */
    stats* best_info; /* info of best iterate               */
	pfloat global_U;
	pfloat global_L;
	
	/* Tmp data */
	char* tmp_bool_node_id;
	pfloat* tmp_int_node_id;
	idxint iter;

	/* Stored pointers to prevent memory leaks */
	pfloat* Gpr_new;
	idxint* Gjc_new;
	idxint* Gir_new;
	pfloat* h_new;

} ecos_bb_pwork;

ecos_bb_pwork* ECOS_BB_setup(
    idxint n, idxint m, idxint p, 
    idxint l, idxint ncones, idxint* q,
    pfloat* Gpr, idxint* Gjc, idxint* Gir,
    pfloat* Apr, idxint* Ajc, idxint* Air,
    pfloat* c, pfloat* h, pfloat* b, 
    idxint num_bool_vars, idxint* bool_vars_idx,
    idxint num_int_vars, idxint* int_vars_idx);

idxint ECOS_BB_solve(ecos_bb_pwork* prob);

void ECOS_BB_cleanup(ecos_bb_pwork* prob, idxint num_vars_keep);

void updateDataEntry_h(ecos_bb_pwork* w, idxint idx, pfloat value);

void updateDataEntry_c(ecos_bb_pwork* w, idxint idx, pfloat value);

/* Calculate the offset into the node_id array */
static inline char* get_bool_node_id(idxint idx, ecos_bb_pwork* prob){
    return &prob->bool_node_ids[prob->num_bool_vars * idx];
}

static inline pfloat* get_int_node_id(idxint idx, ecos_bb_pwork* prob){
    return &prob->int_node_ids[prob->num_int_vars * idx * 2];
}

static inline pfloat abs_2(pfloat number){
	return number < 0.0 ? -number : number;
}

static inline pfloat pfloat_round(pfloat number){
    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

static inline pfloat pfloat_ceil(pfloat number){
	return (pfloat) (number < 0 ? (int) number : (int) (number+(1-MI_INT_TOL)) );
}

static inline pfloat pfloat_floor(pfloat number){
    return (pfloat) (number < 0 ? (int) (number-(1-MI_INT_TOL)) : (int) number);
}

static inline idxint float_eqls(pfloat a, pfloat b){
	return abs_2(a - b) < MI_INT_TOL;
}

#ifdef __cplusplus
}
#endif

#endif
