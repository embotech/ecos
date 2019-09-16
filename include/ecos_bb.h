/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2012-2015 A. Domahidi [domahidi@embotech.com],
 * Automatic Control Lab, ETH Zurich & embotech GmbH, Zurich, Switzerland.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
  * The branch and bound module is (c) Han Wang, Stanford University,
  * [hanwang2@stanford.edu]
  * 
  * Extended with improved branching rules by Pascal Lüscher, student of FHNW
  * [luescherpascal@gmail.com]
  */

#ifndef __ecos_bb_H__
#define __ecos_bb_H__

#include "ecos.h"
#include "spla.h"
#include "glblopts.h"

/* Print verbosity */
#define MI_PRINTLEVEL (1)

/* ecos_bb configuration settings */
#define MI_ABS_EPS (1E-6)
#define MI_REL_EPS (1E-3)
#define MI_MAXITER (1000)
#define MI_INT_TOL (FTOL_INACC)

/* Flags */
#define MI_SOLVED_NON_BRANCHABLE (3)
#define MI_SOLVED_BRANCHABLE (2)
#define MI_NOT_SOLVED (1)
#define MI_FREE (0)

#define MI_ONE (1)
#define MI_ZERO (0)
#define MI_STAR (-1)

/*** Exit flags ***/
/*ECOS_BB found optimal solution*/
#define MI_OPTIMAL_SOLN (ECOS_OPTIMAL)
/*ECOS_BB proved problem is infeasible*/
#define MI_INFEASIBLE (ECOS_PINF)
/*ECOS_BB proved problem is unbounded*/
#define MI_UNBOUNDED (ECOS_DINF)
/*ECOS_BB hit maximum iterations but a feasible solution was found and the best seen feasible solution was returned*/
#define MI_MAXITER_FEASIBLE_SOLN (ECOS_OPTIMAL + ECOS_INACC_OFFSET)
/*ECOS_BB hit maximum iterations without finding a feasible solution*/
#define MI_MAXITER_NO_SOLN (ECOS_PINF + ECOS_INACC_OFFSET)
/*ECOS_BB hit maximum iterations without finding a feasible solution that was unbounded*/
#define MI_MAXITER_UNBOUNDED (ECOS_DINF + ECOS_INACC_OFFSET)

/* Max integer and all smaller integer representable by single precision */
#define MAX_FLOAT_INT (8388608)

#ifdef __cplusplus
extern "C"
{
#endif

	enum BRANCHING_STRATEGY
	{
		BRANCHING_STRATEGY_MOST_INFEASIBLE = 0,
		BRANCHING_STRATEGY_STRONG_BRANCHING = 1,
		BRANCHING_STRATEGY_PSEUDOCOST_BRANCHING = 2,
		BRANCHING_STRATEGY_RELIABILITY = 3,
		BRANCHING_STRATEGY_RANDOM = 4
	};

	enum NODE_SELECTION_METHOD
	{
		BREADTH_FIRST = 0,
		DIVE_LOWER_NODE = 1,
		DIVE_UPPER_NODE = 2,
	};

	typedef struct settings_bb
	{
		idxint maxit;		/* maximum number of iterations         */
		idxint verbose;		/* verbosity bool for PRINTLEVEL < 3    */
		pfloat abs_tol_gap; /* termination criteria |U-L|    		*/
		pfloat rel_tol_gap; /* termination criteria for |U-L|/|L| < 3    */
		pfloat integer_tol; /* integer rounding tolerance */
		enum BRANCHING_STRATEGY branching_strategy;
		idxint reliable_eta; /* number of pseudocost values needed for costs to be reliable */
		enum NODE_SELECTION_METHOD node_selection_method;
	} settings_bb;

	typedef struct node
	{
		char status;
		pfloat L;
		pfloat U;
		pfloat relaxation;
		idxint split_idx;
		pfloat split_val;
		idxint prev_split_idx;
		pfloat prev_split_val;
		pfloat prev_relaxation;
		int up_branch_node;
	} node;

	/* Wrapper for mixed integer module */
	typedef struct ecos_bb_pwork
	{
		/* Mixed integer data */
		idxint num_bool_vars;
		idxint num_int_vars;

		node *nodes;
		char *bool_node_ids;
		pfloat *int_node_ids;

		idxint *bool_vars_idx;
		idxint *int_vars_idx;

		/* ECOS data */
		pwork *ecos_prob;

		/* Modified pointers to ecos internals */
		/* Use these to edit or reset the h variables */
		spmat *A;
		spmat *G;
		pfloat *c;
		pfloat *b;
		pfloat *h;

		/* best iterate seen so far */
		/* variables */
		pfloat *x;   /* primal variables                    */
		pfloat *y;   /* multipliers for equality constaints */
		pfloat *z;   /* multipliers for conic inequalities  */
		pfloat *s;   /* slacks for conic inequalities       */
		pfloat kap;  /* kappa (homogeneous embedding)       */
		pfloat tau;  /* tau (homogeneous embedding)         */
		stats *info; /* info of best iterate               */
		pfloat global_U;
		pfloat global_L;

		/* Tmp data */
		char *tmp_bool_node_id;
		pfloat *tmp_int_node_id;
		idxint iter;
		idxint dive_node_id;

		/* Tmp nodes used for strong branching */
		char *tmp_branching_bool_node_id;
		pfloat *tmp_branching_int_node_id;

		/* Pseudocost branching values */
		pfloat *pseudocost_bin_sum;
		pfloat *pseudocost_int_sum;
		idxint *pseudocost_bin_cnt;
		idxint *pseudocost_int_cnt;

		/* Stored pointers to prevent memory leaks */
		pfloat *Gpr_new;
		idxint *Gjc_new;
		idxint *Gir_new;
		pfloat *h_new;

		/* settings struct */
		settings *ecos_stgs;
		settings_bb *stgs;
		idxint default_settings;

	} ecos_bb_pwork;

	ecos_bb_pwork *ECOS_BB_setup(
		idxint n, idxint m, idxint p,
		idxint l, idxint ncones, idxint *q, idxint nex,
		pfloat *Gpr, idxint *Gjc, idxint *Gir,
		pfloat *Apr, idxint *Ajc, idxint *Air,
		pfloat *c, pfloat *h, pfloat *b,
		idxint num_bool_vars, idxint *bool_vars_idx,
		idxint num_int_vars, idxint *int_vars_idx,
		settings_bb *stgs);

	idxint ECOS_BB_solve(ecos_bb_pwork *prob);

	void ECOS_BB_cleanup(ecos_bb_pwork *prob, idxint num_vars_keep);

	void updateDataEntry_h(ecos_bb_pwork *w, idxint idx, pfloat value);

	void updateDataEntry_c(ecos_bb_pwork *w, idxint idx, pfloat value);

	settings_bb *get_default_ECOS_BB_settings();

	/* Calculate the offset into the node_id array */
	static char *get_bool_node_id(idxint idx, ecos_bb_pwork *prob)
	{
		return &prob->bool_node_ids[prob->num_bool_vars * idx];
	}

	static pfloat *get_int_node_id(idxint idx, ecos_bb_pwork *prob)
	{
		return &prob->int_node_ids[prob->num_int_vars * idx * 2];
	}

	static pfloat abs_2(pfloat number)
	{
		return number < 0.0 ? -number : number;
	}

	static pfloat pfloat_round(pfloat number)
	{
		return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
	}

	static pfloat pfloat_ceil(pfloat number, pfloat integer_tol)
	{
		return (pfloat)(number < 0 ? (int)number : (int)(number + (1 - integer_tol)));
	}

	static pfloat pfloat_floor(pfloat number, pfloat integer_tol)
	{
		return (pfloat)(number < 0 ? (int)(number - (1 - integer_tol)) : (int)number);
	}

	static idxint float_eqls(pfloat a, pfloat b, pfloat integer_tol)
	{
		return abs_2(a - b) < integer_tol;
	}

#ifdef __cplusplus
}
#endif

#endif
