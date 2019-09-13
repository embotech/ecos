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

#include "glblopts.h"
#include "ecos.h"
#include "ecos_bb.h"
#include "math.h"
#include "stdlib.h"
#include "splamm.h"

/* CHOOSE RIGHT MEMORY MANAGER ----------------------------------------- */
#ifdef MATLAB_MEX_FILE
#define MALLOC mxMalloc
#define FREE mxFree
#define CALLOC mxCalloc
#else
#define MALLOC malloc
#define FREE free
#define CALLOC calloc
#endif

static int contains(idxint idx, idxint num_int, idxint *bool_vars_idx)
{
    idxint i;
    idxint ans = 0;
    for (i = 0; i < num_int; ++i)
    {
        ans += (bool_vars_idx[i] == idx);
    }
    return (int)ans;
}

/* Augments the G and b arrays to take lb and ub constraints
 * for all of the variables marked integer
 * USES MALLOC
 */
static void socp_to_ecos_bb(
    idxint num_bool_vars, idxint *bool_vars_idx,
    idxint num_int_vars, idxint *int_vars_idx,
    idxint n, idxint m,
    pfloat *Gpr_in, idxint *Gjc_in, idxint *Gir_in,
    pfloat *Gpr_out, idxint *Gjc_out, idxint *Gir_out,
    pfloat *h_in, pfloat *h_out)
{
    idxint i, j, k;

    /* First map in the column pointers, remember Gir_out[0]=0*/
    for (i = 0; i <= n; ++i)
    {
        Gjc_out[i] = Gjc_in[i];
    }

    j = 0;
    for (i = 0; i < n; ++i)
    {
        if (contains(i, num_bool_vars, bool_vars_idx))
        {
            Gpr_out[Gjc_out[i]] = -1;
            Gpr_out[Gjc_out[i] + 1] = 1;

            Gir_out[Gjc_out[i]] = 2 * j;
            Gir_out[Gjc_out[i] + 1] = 2 * j + 1;

            /* Set lower bound to 0*/
            h_out[2 * j] = 0;

            /* Set upper bound to 1*/
            h_out[2 * j + 1] = 1;

            /* Expand the following arrays */
            for (k = (i + 1); k <= n; ++k)
            {
                Gjc_out[k] += 2;
            }

            /* Now fill in the array */
            for (k = 0; k < (Gjc_in[i + 1] - Gjc_in[i]); ++k)
            {
                Gpr_out[Gjc_out[i] + 2 + k] = Gpr_in[Gjc_in[i] + k];
                Gir_out[Gjc_out[i] + 2 + k] = Gir_in[Gjc_in[i] + k] + 2 * num_bool_vars + 2 * num_int_vars;
            }

            ++j;
        }
        else if (contains(i, num_int_vars, int_vars_idx))
        {
            Gpr_out[Gjc_out[i]] = -1;
            Gpr_out[Gjc_out[i] + 1] = 1;

            Gir_out[Gjc_out[i]] = 2 * j;
            Gir_out[Gjc_out[i] + 1] = 2 * j + 1;

            /* Set lower bound to 0*/
            h_out[2 * j] = MAX_FLOAT_INT;

            /* Set upper bound to 1*/
            h_out[2 * j + 1] = MAX_FLOAT_INT;

            /* Expand the following arrays */
            for (k = (i + 1); k <= n; ++k)
            {
                Gjc_out[k] += 2;
            }

            /* Now fill in the array */
            for (k = 0; k < (Gjc_in[i + 1] - Gjc_in[i]); ++k)
            {
                Gpr_out[Gjc_out[i] + 2 + k] = Gpr_in[Gjc_in[i] + k];
                Gir_out[Gjc_out[i] + 2 + k] = Gir_in[Gjc_in[i] + k] + 2 * num_bool_vars + 2 * num_int_vars;
            }

            ++j;
        }
        else
        {
            for (k = 0; k < (Gjc_in[i + 1] - Gjc_in[i]); ++k)
            {
                Gpr_out[Gjc_out[i] + k] = Gpr_in[Gjc_in[i] + k];
                Gir_out[Gjc_out[i] + k] = Gir_in[Gjc_in[i] + k] + 2 * num_bool_vars + 2 * num_int_vars;
            }
        }
    }

    /* Copy the remaining entries of b*/
    for (i = 0; i < m; ++i)
    {
        h_out[2 * (num_bool_vars + num_int_vars) + i] = h_in[i];
    }
}

ecos_bb_pwork *ECOS_BB_setup(
    idxint n, idxint m, idxint p,
    idxint l, idxint ncones, idxint *q, idxint nex,
    pfloat *Gpr, idxint *Gjc, idxint *Gir,
    pfloat *Apr, idxint *Ajc, idxint *Air,
    pfloat *c, pfloat *h, pfloat *b,
    idxint num_bool_vars, idxint *bool_vars_idx,
    idxint num_int_vars, idxint *int_vars_idx,
    settings_bb *stgs)
{
    idxint new_G_size;
    ecos_bb_pwork *prob;

#if MI_PRINTLEVEL > 2
    int i;
    PRINTTEXT("\n");
    PRINTTEXT("  ****************************************************************\n");
    PRINTTEXT("  * ECOS_BB: Embedded Conic Solver Branch and bound module       *\n");
    PRINTTEXT("  *                                                              *\n");
    PRINTTEXT("  * NOTE: This module is an extension of ECOS by Domahidi et al. *\n");
    PRINTTEXT("  *                                                              *\n");
    PRINTTEXT("  * (C) Han Wang, Stanford University 2014.                      *\n");
    PRINTTEXT("  *                     Email: hanwang2@stanford.edu             *\n");
    PRINTTEXT("  ****************************************************************\n");
    PRINTTEXT("\n\n");
    PRINTTEXT("PROBLEM SUMMARY:\n");
    PRINTTEXT("   Boolean variables (num_bool_vars): %d\n", (int)num_bool_vars);
    PRINTTEXT("   Integer variables ( num_int_vars): %d\n", (int)num_int_vars);
    PRINTTEXT("- - - - - - - - - - - - - - -\n");
    PRINTTEXT("   Boolean var indices: ");
    for (i = 0; i < num_bool_vars; ++i)
        PRINTTEXT("%u ", (unsigned int)bool_vars_idx[i]);
    PRINTTEXT("\n");
    PRINTTEXT("   Integer var indices: ");
    for (i = 0; i < num_int_vars; ++i)
        PRINTTEXT("%u ", (unsigned int)int_vars_idx[i]);
    PRINTTEXT("\n");
    PRINTTEXT("\n");

    PRINTTEXT("n=%u\n", n);
    PRINTTEXT("m=%u\n", m);
    PRINTTEXT("Gpr=");
    for (i = 0; i < Gjc[n]; i++)
    {
        PRINTTEXT("%f,", Gpr[i]);
    }
    PRINTTEXT("\n");
    PRINTTEXT("Gir=");
    for (i = 0; i < Gjc[n]; i++)
    {
        PRINTTEXT("%u,", Gir[i]);
    }
    PRINTTEXT("\n");
    PRINTTEXT("Gjc=");
    for (i = 0; i <= n; i++)
    {
        PRINTTEXT("%u,", Gjc[i]);
    }
    PRINTTEXT("\n");
#endif

    /* MALLOC the problem's memory*/
    prob = (ecos_bb_pwork *)MALLOC(sizeof(ecos_bb_pwork));

    if (stgs == NULL)
    {
        stgs = get_default_ECOS_BB_settings();
        prob->default_settings = 1;
    }
    else
    {
        prob->default_settings = 0;
    }
    prob->stgs = stgs;

    new_G_size = Gjc[n] + (2 * num_bool_vars) + (2 * num_int_vars);
    prob->Gpr_new = (pfloat *)MALLOC(new_G_size * sizeof(pfloat));
    prob->Gjc_new = (idxint *)MALLOC((n + 1) * sizeof(idxint));
    prob->Gir_new = (idxint *)MALLOC(new_G_size * sizeof(idxint));
    prob->h_new = (pfloat *)MALLOC((m + 2 * num_bool_vars + 2 * num_int_vars) * sizeof(pfloat));

    /* Copy the data and convert it to boolean*/
    socp_to_ecos_bb(num_bool_vars, bool_vars_idx,
                    num_int_vars, int_vars_idx,
                    n, m,
                    Gpr, Gjc, Gir,
                    prob->Gpr_new, prob->Gjc_new, prob->Gir_new,
                    h, prob->h_new);
    m += 2 * (num_bool_vars + num_int_vars);
    l += 2 * (num_bool_vars + num_int_vars);

    /* MALLOC the initial node's book keeping data #(2*maxIter)*/
    prob->nodes = (node *)CALLOC(prob->stgs->maxit, sizeof(node));

    /* MALLOC the id arrays*/
    prob->bool_node_ids = (char *)MALLOC(prob->stgs->maxit * num_bool_vars * sizeof(char));
    prob->int_node_ids = (pfloat *)MALLOC(2 * prob->stgs->maxit * num_int_vars * sizeof(pfloat));

    /* MALLOC the tmp node id*/
    prob->tmp_bool_node_id = (char *)MALLOC(num_bool_vars * sizeof(char));
    prob->tmp_int_node_id = (pfloat *)MALLOC(2 * num_int_vars * sizeof(pfloat));

    /* MALLOC the tmp nodes used in branching strategies*/
    prob->tmp_branching_bool_node_id = MALLOC(num_bool_vars * sizeof(char));
    prob->tmp_branching_int_node_id = MALLOC(2 * num_int_vars * sizeof(pfloat));

    /* CALLOC pseudocost branching values*/
    prob->pseudocost_bin_cnt = CALLOC(2 * num_bool_vars, sizeof(idxint));
    prob->pseudocost_bin_sum = CALLOC(2 * num_bool_vars, sizeof(pfloat));
    prob->pseudocost_int_cnt = CALLOC(2 * num_int_vars, sizeof(idxint));
    prob->pseudocost_int_sum = CALLOC(2 * num_int_vars, sizeof(pfloat));

    /* Store the pointer to the boolean idx*/
    prob->bool_vars_idx = bool_vars_idx;
    prob->int_vars_idx = int_vars_idx;

    /* MALLOC the best optimal solution's memory*/
    prob->x = (pfloat *)MALLOC(n * sizeof(pfloat));
    prob->y = (pfloat *)MALLOC(p * sizeof(pfloat));
    prob->z = (pfloat *)MALLOC(m * sizeof(pfloat));
    prob->s = (pfloat *)MALLOC(m * sizeof(pfloat));
    prob->info = (stats *)MALLOC(sizeof(stats));

    /* Setup the ecos solver*/
    prob->ecos_prob = ECOS_setup(
        n, m, p, l, ncones, q, nex,
        prob->Gpr_new, prob->Gjc_new, prob->Gir_new,
        Apr, Ajc, Air,
        c, prob->h_new, b);

    /* Store the number of integer variables in the problem*/
    prob->num_bool_vars = num_bool_vars;
    prob->num_int_vars = num_int_vars;

    prob->global_U = ECOS_INFINITY;

    /* offset the h pointer for the user*/
    prob->h = &prob->ecos_prob->h[2 * (num_bool_vars + num_int_vars)];

    /* Map the other variables*/
    prob->A = prob->ecos_prob->A;
    prob->G = prob->ecos_prob->G;
    prob->c = prob->ecos_prob->c;
    prob->b = prob->ecos_prob->b;

    /* switch off ecos prints */
    prob->ecos_prob->stgs->verbose = 0;

    /* settings */
    prob->ecos_stgs = prob->ecos_prob->stgs;

#if MI_PRINTLEVEL > 2

#if PRINTLEVEL > 2
    PRINTTEXT("ECOS G:\n");
    printSparseMatrix(prob->ecos_prob->G);
#endif

    PRINTTEXT("ECOS h: ");
    for (i = 0; i < m; ++i)
    {
        PRINTTEXT("%.2f ", prob->ecos_prob->h[i]);
    }
    PRINTTEXT("\n");

    PRINTTEXT("\n");
#endif

    return prob;
}

/* Performs the same function as ecos_cleanup*/
void ECOS_BB_cleanup(ecos_bb_pwork *prob, idxint num_vars_keep)
{
    /* FREE solver memory*/
    ECOS_cleanup(prob->ecos_prob, num_vars_keep);
    FREE(prob->tmp_bool_node_id);
    FREE(prob->tmp_int_node_id);
    FREE(prob->tmp_branching_bool_node_id);
    FREE(prob->tmp_branching_int_node_id);
    FREE(prob->Gpr_new);
    FREE(prob->Gjc_new);
    FREE(prob->Gir_new);
    FREE(prob->h_new);
    FREE(prob->nodes);
    FREE(prob->bool_node_ids);
    FREE(prob->int_node_ids);
    FREE(prob->x);
    FREE(prob->y);
    FREE(prob->z);
    FREE(prob->s);
    FREE(prob->info);
    FREE(prob->pseudocost_bin_cnt);
    FREE(prob->pseudocost_bin_sum);
    FREE(prob->pseudocost_int_cnt);
    FREE(prob->pseudocost_int_sum);
    if (prob->default_settings)
    {
        FREE(prob->stgs);
    }
    FREE(prob);
}

settings_bb *get_default_ECOS_BB_settings()
{
    settings_bb *stgs = (settings_bb *)MALLOC(sizeof(settings_bb));
    stgs->maxit = MI_MAXITER;
    stgs->verbose = 1;
    stgs->abs_tol_gap = MI_ABS_EPS;
    stgs->rel_tol_gap = MI_REL_EPS;
    stgs->integer_tol = MI_INT_TOL;
    stgs->branching_strategy = BRANCHING_STRATEGY_RELIABILITY;
    stgs->reliable_eta = 6;
    stgs->node_selection_method = DIVE_LOWER_NODE;
    return stgs;
};
