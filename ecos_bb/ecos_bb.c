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
#include "equil.h"
#include "spla.h"

#include <string.h>

#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

/* Print utility functions*/
#if MI_PRINTLEVEL > 0
static void print_progress(ecos_bb_pwork *prob)
{
    PRINTTEXT("%u \t%.2f \t\t%.2f \t\t%.2f\n", (int)prob->iter, prob->global_L, prob->global_U, prob->global_U - prob->global_L);
}

static void print_ecos_solution(ecos_bb_pwork *prob)
{
    int i;
    PRINTTEXT("ecos->x: ");
    for (i = 0; i < prob->ecos_prob->n; ++i)
        PRINTTEXT("%.2f ", prob->ecos_prob->x[i]);
    PRINTTEXT("\n");
}

static void print_ecos_xequil(ecos_bb_pwork *prob)
{
#if EQUILIBRATE > 0
    int i;
    PRINTTEXT("ecos->xequil: ");
    for (i = 0; i < prob->ecos_prob->n; ++i)
        PRINTTEXT("%.2f ", prob->ecos_prob->xequil[i]);
#else
    PRINTTEXT("ecos->xequil: 1, (equilibration dissabled) ");
#endif
    PRINTTEXT("\n");
}

static void print_ecos_h(ecos_bb_pwork *prob)
{
    int i;
    PRINTTEXT("ecos->h: ");
    for (i = 0; i < prob->ecos_prob->m; ++i)
        PRINTTEXT("%.2f ", prob->ecos_prob->h[i]);
    PRINTTEXT("\n");
}

static void print_ecos_c(ecos_bb_pwork *prob)
{
    int i;
    PRINTTEXT("ecos->c: ");
    for (i = 0; i < prob->ecos_prob->n; ++i)
        PRINTTEXT("%.2f ", prob->ecos_prob->c[i]);
    PRINTTEXT("\n");
}

static void print_node(ecos_bb_pwork *prob, idxint i)
{
    if (i == -1)
    {
        int j;
        PRINTTEXT("Node info: TMP, Partial id:");
        for (j = 0; j < prob->num_bool_vars; ++j)
            PRINTTEXT("%i ", prob->tmp_bool_node_id[j]);
        PRINTTEXT(" | ");
        for (j = 0; j < prob->num_int_vars; ++j)
            PRINTTEXT("(%.2f, %.2f) ", -prob->tmp_int_node_id[2 * j], prob->tmp_int_node_id[2 * j + 1]);
        PRINTTEXT("\n");
    }
    else
    {
        int j;
        PRINTTEXT("Node info %li: %u : %.2f : %.2f : %u : %.2f Partial id:", (long)i, prob->nodes[i].status, prob->nodes[i].L, prob->nodes[i].U, (int)prob->nodes[i].split_idx, prob->nodes[i].split_val);
        for (j = 0; j < prob->num_bool_vars; ++j)
            PRINTTEXT("%i ", get_bool_node_id(i, prob)[j]);
        PRINTTEXT(" | ");
        for (j = 0; j < prob->num_int_vars; ++j)
            PRINTTEXT("(%.2f, %.2f) ", -get_int_node_id(i, prob)[2 * j], get_int_node_id(i, prob)[2 * j + 1]);
        PRINTTEXT("\n");
    }
}

static void print_stats(ecos_bb_pwork *prob)
{
    PRINTTEXT("\tPcost: %.2f \n", prob->info->pcost);
    PRINTTEXT("\tDcost: %.2f \n", prob->info->dcost);
    PRINTTEXT("\tE_Pcost: %.2f \n", prob->ecos_prob->info->pcost);
    PRINTTEXT("\tE_Dcost: %.2f \n", prob->ecos_prob->info->dcost);
}
#endif

static void branch(idxint curr_node_idx, ecos_bb_pwork *prob)
{
    idxint i, split_idx = prob->nodes[curr_node_idx].split_idx;

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose)
    {
        PRINTTEXT("Branching->\t");
        print_node(prob, curr_node_idx);
    }
#endif

    /* Create right node*/
    prob->nodes[prob->iter].L = prob->nodes[curr_node_idx].L;
    prob->nodes[prob->iter].U = prob->nodes[curr_node_idx].U;
    prob->nodes[prob->iter].status = MI_NOT_SOLVED;
    prob->nodes[prob->iter].prev_split_idx = split_idx;
    prob->nodes[prob->iter].prev_split_val = prob->nodes[curr_node_idx].split_val;
    prob->nodes[prob->iter].prev_relaxation = prob->nodes[curr_node_idx].relaxation;
    prob->nodes[prob->iter].up_branch_node = 1;

    prob->nodes[curr_node_idx].prev_split_idx = split_idx;
    prob->nodes[curr_node_idx].prev_split_val = prob->nodes[curr_node_idx].split_val;
    prob->nodes[curr_node_idx].prev_relaxation = prob->nodes[curr_node_idx].relaxation;
    prob->nodes[curr_node_idx].up_branch_node = 0;

    /* Copy over the node id*/
    for (i = 0; i < prob->num_bool_vars; ++i)
        get_bool_node_id(prob->iter, prob)[i] = get_bool_node_id(curr_node_idx, prob)[i];
    for (i = 0; i < prob->num_int_vars * 2; ++i)
        get_int_node_id(prob->iter, prob)[i] = get_int_node_id(curr_node_idx, prob)[i];

    if (split_idx < prob->num_bool_vars)
    {
        get_bool_node_id(curr_node_idx, prob)[split_idx] = MI_ZERO;
        get_bool_node_id(prob->iter, prob)[split_idx] = MI_ONE;
    }
    else
    {
        split_idx -= prob->num_bool_vars;

        /* Left branch constrain UB */
        get_int_node_id(curr_node_idx, prob)[split_idx * 2 + 1] =
            pfloat_floor(prob->nodes[curr_node_idx].split_val, prob->stgs->integer_tol);

        /* Right branch constrain LB */
        get_int_node_id(prob->iter, prob)[split_idx * 2] =
            -pfloat_ceil(prob->nodes[curr_node_idx].split_val, prob->stgs->integer_tol);
    }

    prob->nodes[curr_node_idx].status = MI_NOT_SOLVED;

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose)
    {
        PRINTTEXT(" Left-> \t ");
        print_node(prob, curr_node_idx);
        PRINTTEXT(" Right->\t ");
        print_node(prob, prob->iter);
    }
#endif
}

/*
 * Function to return the next node to be expanded
 */
static idxint get_next_node(ecos_bb_pwork *prob)
{
    if (prob->stgs->node_selection_method == DIVE_LOWER_NODE && prob->nodes[prob->dive_node_id].status == MI_SOLVED_BRANCHABLE)
    {
        return prob->dive_node_id;
    }

    if (prob->stgs->node_selection_method == DIVE_UPPER_NODE && prob->nodes[prob->iter].status == MI_SOLVED_BRANCHABLE)
    {
        return prob->iter;
    }

    idxint i;
    idxint next_node = -1;
    pfloat L = ECOS_INFINITY;
    for (i = 0; i <= prob->iter; ++i)
    {
        if (prob->nodes[i].status == MI_SOLVED_BRANCHABLE && prob->nodes[i].L < L && prob->nodes[i].L < prob->global_U)
        {
            next_node = i;
            L = prob->nodes[i].L;
        }
    }
    prob->dive_node_id = next_node;
    return next_node;
}

static pfloat get_global_L(ecos_bb_pwork *prob)
{
    idxint i;
    pfloat L = ECOS_INFINITY;
    for (i = 0; i <= prob->iter; ++i)
        L = MIN(L, prob->nodes[i].L);
    return L;
}

/*
 * Function to return the next var to split on
 * using most infeasible branching
 */
static void get_branch_var_most_infeasible(ecos_bb_pwork *prob, idxint *split_idx, pfloat *split_val)
{
    idxint i;
    pfloat x, y, d, ambiguity;
    d = 1.0;
    for (i = 0; i < (prob->num_bool_vars + prob->num_int_vars); ++i)
    {
        if (i < prob->num_bool_vars)
        {
            y = prob->ecos_prob->x[prob->bool_vars_idx[i]];
            x = y;
        }
        else
        {
            y = prob->ecos_prob->x[prob->int_vars_idx[i - prob->num_bool_vars]];
            x = y - pfloat_floor(y, prob->stgs->integer_tol);
        }
        ambiguity = abs_2(x - 0.5);
        if (ambiguity < d)
        {
            *split_idx = i;
            *split_val = y;
            d = ambiguity;
        }
    }
#if MI_PRINTLEVEL > 1
    PRINTTEXT("split_idx:%u, split_val:%f\n", *split_idx, *split_val);
#endif
}

/*
 * Function to return the next var to split on
 * using random branching
 */
static void get_branch_var_random(ecos_bb_pwork *prob, idxint *split_idx, pfloat *split_val)
{
    while (1)
    {
        idxint i = rand() % (prob->num_bool_vars + prob->num_int_vars);
        *split_val = i < prob->num_bool_vars
                         ? prob->ecos_prob->x[prob->bool_vars_idx[i]]
                         : prob->ecos_prob->x[prob->int_vars_idx[i - prob->num_bool_vars]];
        *split_idx = i;
        if (float_eqls(*split_val, pfloat_floor(*split_val, prob->stgs->integer_tol), prob->stgs->integer_tol) ||
            float_eqls(*split_val, pfloat_ceil(*split_val, prob->stgs->integer_tol), prob->stgs->integer_tol))
        {
            continue;
        }

        break;
    }
#if MI_PRINTLEVEL > 1
    PRINTTEXT("split_idx:%u, split_val:%f\n", *split_idx, *split_val);
#endif
}

/*
 * Updates the solver's lb and ub contraints for integer variables
 * to the bounds specified by the node
 */
static void set_prob(ecos_bb_pwork *prob, char *bool_node_id, pfloat *int_node_id)
{
    idxint i;
    for (i = 0; i < prob->num_bool_vars; ++i)
    {
        switch (bool_node_id[i])
        {
        case MI_ONE:
            ecos_updateDataEntry_h(prob->ecos_prob, 2 * i, -1.0);    /*ecos_prob->h[2*i] = -(1.0); // -x <= -1 <=> x >= 1*/
            ecos_updateDataEntry_h(prob->ecos_prob, 2 * i + 1, 1.0); /*ecos_prob->h[2*i + 1] = 1.0;*/
            break;
        case MI_ZERO:
            ecos_updateDataEntry_h(prob->ecos_prob, 2 * i, 0.0);     /*ecos_prob->h[2*i] = 0.0;*/
            ecos_updateDataEntry_h(prob->ecos_prob, 2 * i + 1, 0.0); /*ecos_prob->h[2*i + 1] = 0.0;*/
            break;
        case MI_STAR:
            ecos_updateDataEntry_h(prob->ecos_prob, 2 * i, 0.0);     /*ecos_prob->h[2*i] = 0.0;*/
            ecos_updateDataEntry_h(prob->ecos_prob, 2 * i + 1, 1.0); /*ecos_prob->h[2*i + 1] = 1.0;*/
            break;
        default:
#if MI_PRINTLEVEL > 1
            PRINTTEXT("Illegal boolean setting arguments passed: %u \n", bool_node_id[i]);
#endif
            break;
        }
    }

    /* Set integer bounds */
    for (i = 0; i < prob->num_int_vars; ++i)
    {
        ecos_updateDataEntry_h(prob->ecos_prob, 2 * (i + prob->num_bool_vars), int_node_id[2 * i]);
        ecos_updateDataEntry_h(prob->ecos_prob, 2 * (i + prob->num_bool_vars) + 1, int_node_id[2 * i + 1]);
    }

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose)
    {
        print_ecos_h(prob);
    }
#endif
}

static char is_infeasible(idxint ecos_result)
{
    if (ecos_result == ECOS_DINF ||
        ecos_result == ECOS_PINF ||
        ecos_result == ECOS_DINF + ECOS_INACC_OFFSET ||
        ecos_result == ECOS_PINF + ECOS_INACC_OFFSET)
    {
        return 1;
    }

    return 0;
}

/*
 * returns the score for down / up branching.
 * Using the product, same method is used in SCIP
 */
static pfloat get_score(const pfloat delta_q_down, const pfloat delta_q_up)
{

    const pfloat epsilon = 0.00000001;
    return (delta_q_down > epsilon ? delta_q_down : epsilon) * (delta_q_up > epsilon ? delta_q_up : epsilon);

    // alternative score function using weighted avg
    // const pfloat mu = 1.0 / 6.0;
    // return (1 - mu) * min(delta_q_down, delta_q_up) + mu * max(delta_q_down, delta_q_up);
}

/*
 * Creates a new array with length n with random sorted idxints
 */
static idxint *get_shuffled_array(const idxint n)
{
    idxint *ret_val = MALLOC(n * sizeof(idxint));
    idxint i;
    for (i = 0; i < n; ++i)
    {
        ret_val[i] = i;
    }

    for (i = 0; i < n - 1; i++)
    {
        idxint j = i + rand() / (RAND_MAX / (n - i) + 1);
        const idxint t = ret_val[j];
        ret_val[j] = ret_val[i];
        ret_val[i] = t;
    }
    return ret_val;
}

/**
 * Initializes strong branching by copying current state in tmp_branching nodes
 */
static void initialize_strong_branching(ecos_bb_pwork *problem, idxint node_idx)
{
    // copy away tmp bool and int node since it's needed by another function -.-
    const int bool_node_size = sizeof(char) * problem->num_bool_vars;
    const int int_node_size = sizeof(pfloat) * 2 * problem->num_int_vars;

    // copy over current state to tmp_branching_nodes

    memcpy(problem->tmp_branching_bool_node_id, get_bool_node_id(node_idx, problem), bool_node_size);
    memcpy(problem->tmp_branching_int_node_id, get_int_node_id(node_idx, problem), int_node_size);
}

/*
 * Calculates the relaxation (and ecos_result) for the tmp_branching node
 */
static void calc_tmp_branching_problem(ecos_bb_pwork *problem, pfloat *relaxation, idxint *ecos_result)
{
    set_prob(problem, problem->tmp_branching_bool_node_id, problem->tmp_branching_int_node_id);
    *ecos_result = ECOS_solve(problem->ecos_prob);
    *relaxation = eddot(problem->ecos_prob->n, problem->ecos_prob->x, problem->ecos_prob->c);
}

/*
 * Checks if the ecos result is infeasible, if so it fixes the variable to the other_val
 */
static int check_infeasible_bool_var(ecos_bb_pwork *problem, const idxint ecos_result, const pfloat relaxation, const idxint node_idx,
                                     const char other_val, const idxint var_idx, idxint *split_idx, pfloat *split_val,
                                     const pfloat current_value)
{
    if (is_infeasible(ecos_result) || relaxation > problem->global_U)
    {
        //fix value to MI_ONE, no branching needed
        get_bool_node_id(node_idx, problem)[var_idx] = other_val;
        problem->tmp_branching_bool_node_id[var_idx] = other_val;
#if PRINTLEVEL >= 3
        PRINTTEXT("Found infeasible node path for down-branching in node [%d] in var [%d]\n", node_idx, i);
#endif
        // prevent to have no branching variable set because of continuing
        if (*split_idx == -1)
        {
            *split_idx = var_idx;
            *split_val = current_value;
        }
        return 1;
    }
    return 0;
}

/*
 * Checks if the ecos result is infeasible, if so it fixes the variable
 */
static int check_infeasible_int_var(ecos_bb_pwork *problem, const idxint ecos_result, const pfloat relaxation, const idxint node_idx,
                                    const int fix_lb, const idxint var_idx, idxint *split_idx, pfloat *split_val,
                                    const pfloat current_value)
{
    if (is_infeasible(ecos_result) || relaxation > problem->global_U)
    {
        // fix lower bound
        const pfloat fix_val = fix_lb
                                   ? -pfloat_ceil(current_value, problem->stgs->integer_tol)
                                   : pfloat_floor(current_value, problem->stgs->integer_tol);
        get_int_node_id(node_idx, problem)[2 * var_idx + (fix_lb ? 0 : 1)] = fix_val;
        problem->tmp_branching_int_node_id[2 * var_idx + (fix_lb ? 0 : 1)] = fix_val;
#if PRINTLEVEL >= 3
        PRINTTEXT("Found infeasible node path for down-branching in node [%d] in var [%d]\n", node_idx, i);
#endif
        // prevent to have no branching variable set because of continuing
        if (*split_idx == -1)
        {
            *split_idx = var_idx + problem->num_bool_vars;
            *split_val = current_value;
        }
        return 1;
    }
    return 0;
}

/*
 * calculates the q_down and q_up using strong-branching for a bool var
 * returns 1 if one of the paths is infeasible
 */
static int strong_branch_bool_var(ecos_bb_pwork *problem, idxint *split_idx, pfloat *split_val, const idxint node_idx,
                                  pfloat *q_down, pfloat *q_up, const idxint i, const pfloat current_value)
{
    // save val before branching
    const char orig_val = problem->tmp_branching_bool_node_id[i];
    idxint ecos_result;
    // branch down
    problem->tmp_branching_bool_node_id[i] = MI_ZERO;
    calc_tmp_branching_problem(problem, q_down, &ecos_result);
    if (check_infeasible_bool_var(problem, ecos_result, *q_down, node_idx, MI_ONE, i, split_idx, split_val, current_value))
    {
        return 1;
    }
    // branch up
    problem->tmp_branching_bool_node_id[i] = MI_ONE;
    calc_tmp_branching_problem(problem, q_up, &ecos_result);
    if (check_infeasible_bool_var(problem, ecos_result, *q_up, node_idx, MI_ZERO, i, split_idx, split_val, current_value))
    {
        return 1;
    }

    problem->tmp_branching_bool_node_id[i] = orig_val;
    return 0;
}

/*
 * calculates the q_down and q_up using strong-branching for an int var
 * returns 1 if one of the paths is infeasible
 */
static int strong_branch_int_var(ecos_bb_pwork *problem, idxint *split_idx, pfloat *split_val, idxint node_idx, pfloat *q_down, pfloat *q_up, idxint i, pfloat current_value)
{
    const idxint int_idx = i - problem->num_bool_vars;
    idxint ecos_result;
    // branch down (set upper bound)
    pfloat orig_val = problem->tmp_branching_int_node_id[2 * int_idx + 1];
    problem->tmp_branching_int_node_id[2 * int_idx + 1] = pfloat_floor(current_value, problem->stgs->integer_tol);
    calc_tmp_branching_problem(problem, q_down, &ecos_result);
    problem->tmp_branching_int_node_id[2 * int_idx + 1] = orig_val;
    if (check_infeasible_int_var(problem, ecos_result, *q_down, node_idx, 1, int_idx, split_idx, split_val, current_value))
    {
        return 1;
    }

    // branch up (set lower bound)
    orig_val = problem->tmp_branching_int_node_id[2 * int_idx];
    problem->tmp_branching_int_node_id[2 * int_idx] = -pfloat_ceil(current_value, problem->stgs->integer_tol);
    calc_tmp_branching_problem(problem, q_up, &ecos_result);

    problem->tmp_branching_int_node_id[2 * int_idx] = orig_val;
    if (check_infeasible_int_var(problem, ecos_result, *q_up, node_idx, 0, int_idx, split_idx, split_val, current_value))
    {
        return 1;
    }
    return 0;
}

/*
 * Function to return the next var to split on
 */
static void get_branch_var_strong_branching(ecos_bb_pwork *problem, idxint *split_idx, pfloat *split_val, idxint node_idx)
{
    /* id of either bool index or int index */
    idxint var_id;
    /* best cost based on get_score */
    pfloat best_score = -ECOS_INFINITY;

    initialize_strong_branching(problem, node_idx);

    /* calc current lp relaxation and save x vector */
    idxint ecos_result;
    pfloat q;
    calc_tmp_branching_problem(problem, &q, &ecos_result);
    pfloat *x_values = MALLOC(sizeof(pfloat) * problem->ecos_prob->n);
    memcpy(x_values, problem->ecos_prob->x, sizeof(pfloat) * problem->ecos_prob->n);

    /* the lp-relaxation for the down branching solution */
    pfloat q_down;
    /* the lp-relaxation for the up branching solution */
    pfloat q_up;

    pfloat best_q_down = 0.0, best_q_up = 0.0;

    idxint orig_maxit = problem->ecos_prob->stgs->maxit;
    idxint orig_maxit_bk = problem->ecos_prob->stgs->max_bk_iter;
    /* set split idx to something impossible */
    *split_idx = -1;

    idxint *random_idx = get_shuffled_array(problem->num_bool_vars + problem->num_int_vars);
    idxint rand_i;
    for (rand_i = 0; rand_i < problem->num_bool_vars + problem->num_int_vars; ++rand_i)
    {
        idxint const i = random_idx[rand_i];

        if (i < problem->num_bool_vars)
        {
            var_id = problem->bool_vars_idx[i];
        }
        else
        {

            var_id = problem->int_vars_idx[i - problem->num_bool_vars];
        }
        pfloat current_value = x_values[var_id];
        // if already int, skip candidate
        if (float_eqls(current_value, pfloat_floor(current_value, problem->stgs->integer_tol), problem->stgs->integer_tol) ||
            float_eqls(current_value, pfloat_ceil(current_value, problem->stgs->integer_tol), problem->stgs->integer_tol))
        {
            continue;
        }

        if (i < problem->num_bool_vars)
        {
            if (strong_branch_bool_var(problem, split_idx, split_val, node_idx, &q_down, &q_up, i, current_value))
                continue;
        }
        else
        {
            if (strong_branch_int_var(problem, split_idx, split_val, node_idx, &q_down, &q_up, i, current_value))
                continue;
        }

        const pfloat curr_score = get_score(q_down - q, q_up - q);
        if (curr_score > best_score)
        {
            *split_idx = i;
            *split_val = current_value;
            best_q_down = q_down;
            best_q_up = q_up;

            best_score = curr_score;
        }
    }
    free(x_values);
    free(random_idx);
    problem->ecos_prob->stgs->maxit = orig_maxit;
    problem->ecos_prob->stgs->max_bk_iter = orig_maxit_bk;

#if MI_PRINTLEVEL >= 3
    PRINTTEXT("split_idx: %u, split_val: %f,q: %f, q_down: %f, q_up: %f\n", *split_idx, *split_val, q, best_q_down, best_q_up);
#endif
}

/*
 * returns the average pseudocost values
 */
static pfloat avg_pseudocost(pfloat *pseudocost_sums, idxint *pseudocost_cnts, const int pseudocost_sums_cnt, const char up)
{
    pfloat sum = 0.0;
    idxint cnt = 0;
    idxint i = 0;
    for (i = 0; i < pseudocost_sums_cnt; ++i)
    {
        if (pseudocost_cnts[2 * i + up] > 0)
        {
            sum += pseudocost_sums[2 * i + up];
            cnt += pseudocost_cnts[2 * i + up];
        }
    }
    if (cnt > 0)
    {
        return sum / cnt;
    }
    return 1.0;
}

/*
 * calculates the psi values for pseudocosts
 */
static void set_pseudocost_psi(ecos_bb_pwork *problem, const idxint i, pfloat *down_psi, pfloat *up_psi)
{
    if (i < problem->num_bool_vars)
    {
        if (problem->pseudocost_bin_cnt[2 * i] == 0)
        {
            *down_psi = avg_pseudocost(problem->pseudocost_bin_sum, problem->pseudocost_bin_cnt, problem->num_bool_vars, 0);
        }
        else
        {
            *down_psi = (problem->pseudocost_bin_sum[2 * i] / problem->pseudocost_bin_cnt[2 * i]);
        }

        if (problem->pseudocost_bin_cnt[2 * i + 1] == 0)
        {
            *up_psi = avg_pseudocost(problem->pseudocost_bin_sum, problem->pseudocost_bin_cnt, problem->num_bool_vars, 1);
        }
        else
        {
            *up_psi = (problem->pseudocost_bin_sum[2 * i + 1] / problem->pseudocost_bin_cnt[2 * i + 1]);
        }
    }
    else
    {
        idxint var_idx = i - problem->num_bool_vars;
        if (problem->pseudocost_int_cnt[2 * var_idx] == 0)
        {
            *down_psi = avg_pseudocost(problem->pseudocost_int_sum, problem->pseudocost_int_cnt, problem->num_int_vars, 0);
        }
        else
        {
            *down_psi = (problem->pseudocost_int_sum[2 * var_idx] / problem->pseudocost_int_cnt[2 * var_idx]);
            ;
        }

        if (problem->pseudocost_int_cnt[2 * var_idx + 1] == 0)
        {
            *up_psi = avg_pseudocost(problem->pseudocost_int_sum, problem->pseudocost_int_cnt, problem->num_int_vars, 1);
        }
        else
        {
            *up_psi = (problem->pseudocost_int_sum[2 * var_idx + 1] / problem->pseudocost_int_cnt[2 * var_idx + 1]);
        }
    }
}

/*
 * Function to return the next var to split on
 * using pure pseudocost branching
 */
static void get_branch_var_pseudocost_branching(ecos_bb_pwork *problem, idxint *split_idx, pfloat *split_val, idxint node_idx)
{
    pfloat high_score = -ECOS_INFINITY;
    pfloat up_psi, down_psi;
    idxint var_idx;
    idxint *random_idx = get_shuffled_array(problem->num_bool_vars + problem->num_int_vars);
    idxint rand_i;
    for (rand_i = 0; rand_i < problem->num_bool_vars + problem->num_int_vars; ++rand_i)
    {
        idxint const i = random_idx[rand_i];
        if (i < problem->num_bool_vars)
        {
            var_idx = problem->bool_vars_idx[i];
        }
        else
        {
            var_idx = problem->int_vars_idx[(i - problem->num_bool_vars)];
        }

        pfloat y = problem->ecos_prob->x[var_idx];
        if (float_eqls(y, pfloat_floor(y, problem->stgs->integer_tol), problem->stgs->integer_tol) ||
            float_eqls(y, pfloat_ceil(y, problem->stgs->integer_tol), problem->stgs->integer_tol))
        {
            continue;
        }

        set_pseudocost_psi(problem, i, &down_psi, &up_psi);
        const pfloat down_change_in_var = y - pfloat_floor(y, problem->stgs->integer_tol);
        const pfloat up_change_in_var = 1.0 - down_change_in_var;

        const pfloat score = get_score(down_change_in_var * down_psi, up_change_in_var * up_psi);

        if (score > high_score)
        {
            *split_idx = i;
            *split_val = y;
            high_score = score;
        }
    }
    free(random_idx);
#if MI_PRINTLEVEL > 1
    PRINTTEXT("split_idx:%u, split_val:%f\n", *split_idx, *split_val);
#endif
}

/*
 * Checks if the variable i is reliable enough for pseudocost branching
 */
static int is_reliable(ecos_bb_pwork *problem, const idxint i)
{
    if (i < problem->num_bool_vars)
    {
        return problem->pseudocost_bin_cnt[2 * i] > problem->stgs->reliable_eta &&
               problem->pseudocost_bin_cnt[2 * i + 1] > problem->stgs->reliable_eta;
    }
    return problem->pseudocost_int_cnt[2 * (i - problem->num_bool_vars)] > problem->stgs->reliable_eta &&
           problem->pseudocost_int_cnt[2 * (i - problem->num_bool_vars) + 1] > problem->stgs->reliable_eta;
}

/*
 * Function to return the next var to split on
 * using reliability branching
 */
static void get_branch_var_reliability_branching(ecos_bb_pwork *problem, idxint *split_idx, pfloat *split_val, const idxint node_idx)
{
    pfloat score;
    pfloat high_score = -ECOS_INFINITY;
    pfloat up_psi, down_psi;
    idxint var_idx;

    int strong_initialized = 0;

    idxint ecos_result;
    pfloat q;
    pfloat *x_values = NULL;

    idxint *random_idx = get_shuffled_array(problem->num_bool_vars + problem->num_int_vars);
    idxint rand_i;
    for (rand_i = 0; rand_i < problem->num_bool_vars + problem->num_int_vars; ++rand_i)
    {
        idxint const i = random_idx[rand_i];

        if (i < problem->num_bool_vars)
        {
            var_idx = problem->bool_vars_idx[i];
        }
        else
        {
            var_idx = problem->int_vars_idx[(i - problem->num_bool_vars)];
        }

        const pfloat y = problem->ecos_prob->x[var_idx];
        if (float_eqls(y, pfloat_floor(y, problem->stgs->integer_tol), problem->stgs->integer_tol) ||
            float_eqls(y, pfloat_ceil(y, problem->stgs->integer_tol), problem->stgs->integer_tol))
        {
            continue;
        }

        if (is_reliable(problem, i))
        {
            set_pseudocost_psi(problem, i, &down_psi, &up_psi);

            const pfloat down_change_in_var = y - pfloat_floor(y, problem->stgs->integer_tol);
            const pfloat up_change_in_var = 1.0 - down_change_in_var;

            score = get_score(down_change_in_var * down_psi, up_change_in_var * up_psi);
        }
        else
        {
            if (!strong_initialized)
            {
                initialize_strong_branching(problem, node_idx);
                strong_initialized = 1;

                calc_tmp_branching_problem(problem, &q, &ecos_result);
                x_values = MALLOC(sizeof(pfloat) * problem->ecos_prob->n);
                memcpy(x_values, problem->ecos_prob->x, sizeof(pfloat) * problem->ecos_prob->n);
            }

            var_idx = i;

            pfloat q_down, q_up;
            pfloat *sum_ptr;
            idxint *cnt_ptr;
            pfloat current_value;
            if (i < problem->num_bool_vars)
            {
                current_value = x_values[problem->bool_vars_idx[var_idx]];
                if (strong_branch_bool_var(problem, split_idx, split_val, node_idx, &q_down, &q_up, i, current_value))
                {
                    continue;
                }
                sum_ptr = problem->pseudocost_bin_sum;
                cnt_ptr = problem->pseudocost_bin_cnt;
            }
            else
            {
                var_idx = i - problem->num_bool_vars;
                current_value = x_values[problem->int_vars_idx[var_idx]];
                if (strong_branch_int_var(problem, split_idx, split_val, node_idx, &q_down, &q_up, i, x_values[problem->int_vars_idx[var_idx]]))
                {
                    continue;
                }
                sum_ptr = problem->pseudocost_int_sum;
                cnt_ptr = problem->pseudocost_int_cnt;
            }
            const pfloat var_change_up = pfloat_ceil(current_value, problem->stgs->integer_tol) - current_value;
            const pfloat var_change_down = current_value - pfloat_floor(current_value, problem->stgs->integer_tol);
            sum_ptr[2 * var_idx] += (q_down - q) / var_change_down;
            cnt_ptr[2 * var_idx] += 1;
            sum_ptr[2 * var_idx + 1] += (q_up - q) / var_change_up;
            cnt_ptr[2 * var_idx + 1] += 1;

            score = get_score(q_down - q, q_up - q);
        }

        if (score > high_score)
        {
            *split_idx = i;
            *split_val = y;
            high_score = score;
        }
    }
    free(random_idx);
    free(x_values);
#if MI_PRINTLEVEL > 1
    PRINTTEXT("split_idx:%u, split_val:%f\n", *split_idx, *split_val);
#endif
}

/*
 * Stores the ecos solution to the array inside ecos_bb
 */
static void store_solution(ecos_bb_pwork *prob)
{
    idxint i;
    for (i = 0; i < prob->ecos_prob->n; ++i)
        prob->x[i] = prob->ecos_prob->x[i];
    for (i = 0; i < prob->ecos_prob->p; ++i)
        prob->y[i] = prob->ecos_prob->y[i];
    for (i = 0; i < prob->ecos_prob->m; ++i)
        prob->z[i] = prob->ecos_prob->z[i];
    for (i = 0; i < prob->ecos_prob->m; ++i)
        prob->s[i] = prob->ecos_prob->s[i];
    prob->kap = prob->ecos_prob->kap;
    prob->tau = prob->ecos_prob->tau;
    *(prob->info) = *(prob->ecos_prob->info);
}

/*
 * Loads the ecos_bb solution back into the ecos struct, necessary for the Matlab/Python interface
 */
static void load_solution(ecos_bb_pwork *prob)
{
    idxint i;
    for (i = 0; i < prob->ecos_prob->n; ++i)
        prob->ecos_prob->x[i] = prob->x[i];
    for (i = 0; i < prob->ecos_prob->p; ++i)
        prob->ecos_prob->y[i] = prob->y[i];
    for (i = 0; i < prob->ecos_prob->m; ++i)
        prob->ecos_prob->z[i] = prob->z[i];
    for (i = 0; i < prob->ecos_prob->m; ++i)
        prob->ecos_prob->s[i] = prob->s[i];
    prob->ecos_prob->kap = prob->kap;
    prob->ecos_prob->tau = prob->tau;
    *(prob->ecos_prob->info) = *(prob->info);
}

static void get_bounds(idxint node_idx, ecos_bb_pwork *prob)
{
    idxint i, ret_code, branchable, viable_rounded_sol;
    viable_rounded_sol = 0;

    set_prob(prob, get_bool_node_id(node_idx, prob), get_int_node_id(node_idx, prob));
    ret_code = ECOS_solve(prob->ecos_prob);

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose)
    {
        print_ecos_solution(prob);
    }
    if (ret_code != ECOS_OPTIMAL && ret_code != ECOS_PINF && ret_code != ECOS_DINF)
    {
        PRINTTEXT("1Exit code: %u\n", ret_code);
    }
#endif

    if (ret_code == ECOS_OPTIMAL ||
        ret_code == ECOS_OPTIMAL + ECOS_INACC_OFFSET ||
        ret_code == ECOS_MAXIT ||
        ret_code == ECOS_NUMERICS)
    {
        prob->nodes[node_idx].L = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);

        // if global best known solution is better than our lower
        if (prob->global_U < prob->nodes[node_idx].L)
        {
            prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
            return;
        }

        prob->nodes[node_idx].relaxation = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);

        // update pseudocost_values
        if (prob->nodes[node_idx].prev_split_idx >= 0)
        {
            idxint split_idx = prob->nodes[node_idx].prev_split_idx;
            const idxint offset = prob->nodes[node_idx].up_branch_node ? 1 : 0;

            const pfloat change_in_q = prob->nodes[node_idx].relaxation - prob->nodes[node_idx].prev_relaxation;
            const pfloat change_in_var = prob->nodes[node_idx].up_branch_node
                                             ? pfloat_ceil(prob->nodes[node_idx].prev_split_val, prob->stgs->integer_tol) - prob->nodes[node_idx].prev_split_val
                                             : prob->nodes[node_idx].prev_split_val - pfloat_floor(prob->nodes[node_idx].prev_split_val, prob->stgs->integer_tol);

            pfloat *sum_ptr;
            idxint *cnt_ptr;
            if (split_idx < prob->num_bool_vars)
            {
                sum_ptr = &prob->pseudocost_bin_sum[2 * split_idx + offset];
                cnt_ptr = &prob->pseudocost_bin_cnt[2 * split_idx + offset];
            }
            else
            {
                split_idx -= prob->num_bool_vars;
                sum_ptr = &prob->pseudocost_int_sum[2 * split_idx + offset];
                cnt_ptr = &prob->pseudocost_int_cnt[2 * split_idx + offset];
            }
            *sum_ptr += change_in_q / change_in_var;
            *cnt_ptr += 1;
        }

        /* Figure out if x is already an integer solution
			if the solution had no numerical errors*/
        branchable = 1;
        for (i = 0; i < prob->num_bool_vars; ++i)
        {
            prob->tmp_bool_node_id[i] = (char)pfloat_round(prob->ecos_prob->x[prob->bool_vars_idx[i]]);
            branchable &= float_eqls(prob->ecos_prob->x[prob->bool_vars_idx[i]], (pfloat)prob->tmp_bool_node_id[i], prob->stgs->integer_tol);
        }
        for (i = 0; i < prob->num_int_vars; ++i)
        {
            prob->tmp_int_node_id[2 * i + 1] = pfloat_round(prob->ecos_prob->x[prob->int_vars_idx[i]]);
            prob->tmp_int_node_id[2 * i] = -(prob->tmp_int_node_id[2 * i + 1]);
            branchable &= float_eqls(prob->ecos_prob->x[prob->int_vars_idx[i]], prob->tmp_int_node_id[2 * i + 1], prob->stgs->integer_tol);
        }
        branchable = !branchable;

#if MI_PRINTLEVEL > 1
        if (prob->stgs->verbose)
        {
            PRINTTEXT("Is branchable: %u\n", branchable);
        }
#endif

        if (branchable)
        { /* pfloat_round and check feasibility*/
            switch (prob->stgs->branching_strategy)
            {
            case BRANCHING_STRATEGY_MOST_INFEASIBLE:
                get_branch_var_most_infeasible(prob, &(prob->nodes[node_idx].split_idx), &(prob->nodes[node_idx].split_val));
                break;
            case BRANCHING_STRATEGY_STRONG_BRANCHING:
                get_branch_var_strong_branching(prob, &(prob->nodes[node_idx].split_idx), &(prob->nodes[node_idx].split_val), node_idx);
                break;
            case BRANCHING_STRATEGY_PSEUDOCOST_BRANCHING:
                get_branch_var_pseudocost_branching(prob, &(prob->nodes[node_idx].split_idx), &(prob->nodes[node_idx].split_val), node_idx);
                break;
            case BRANCHING_STRATEGY_RELIABILITY:
                get_branch_var_reliability_branching(prob, &(prob->nodes[node_idx].split_idx), &(prob->nodes[node_idx].split_val), node_idx);
                break;
            case BRANCHING_STRATEGY_RANDOM:
                get_branch_var_random(prob, &(prob->nodes[node_idx].split_idx), &(prob->nodes[node_idx].split_val));
            default:;
            }

            prob->nodes[node_idx].status = MI_SOLVED_BRANCHABLE;

#if MI_PRINTLEVEL > 1
            if (prob->stgs->verbose)
            {
                print_node(prob, -1);
                PRINTTEXT("Rounded Solution:\n");
            }
#endif
            set_prob(prob, prob->tmp_bool_node_id, prob->tmp_int_node_id);
            ret_code = ECOS_solve(prob->ecos_prob);

#if MI_PRINTLEVEL > 1
            if (prob->stgs->verbose)
            {
                print_ecos_solution(prob);
            }
            if (ret_code != ECOS_OPTIMAL && ret_code != ECOS_PINF && ret_code != ECOS_DINF)
            {
                PRINTTEXT("2Exit code: %u\n", ret_code);
            }
#endif
            if (ret_code == ECOS_OPTIMAL)
            {
                /* Use the node's U as tmp storage */
                prob->nodes[node_idx].U = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);
                viable_rounded_sol = 1;
            }
        }
        else
        { /* This is already an integer solution*/
            prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
            prob->nodes[node_idx].U = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);
        }

        if (prob->nodes[node_idx].U < prob->global_U)
        {
#if MI_PRINTLEVEL > 1
            if (prob->stgs->verbose)
            {
                PRINTTEXT("New optimal solution, U: %.2f\n", prob->nodes[node_idx].U);
                print_ecos_xequil(prob);
                print_ecos_c(prob);
                print_ecos_solution(prob);
            }
#endif
            store_solution(prob);
            prob->global_U = prob->nodes[node_idx].U;
        }

        if (viable_rounded_sol)
        {
            /* Reset the node's U back to INF because it was not originally feasible */
            prob->nodes[node_idx].U = ECOS_INFINITY;
        }

        // if global best known solution is better than our lower
        if (prob->global_U < prob->nodes[node_idx].L)
        {
            prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
        }
    }
    else
    { /*Assume node infeasible*/
        prob->nodes[node_idx].L = ECOS_INFINITY;
        prob->nodes[node_idx].U = ECOS_INFINITY;
        prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
    }
}

static idxint should_continue(ecos_bb_pwork *prob, idxint curr_node_idx)
{
    return (prob->global_U - prob->global_L) > prob->stgs->abs_tol_gap && abs_2(prob->global_U / prob->global_L - 1.0) > prob->stgs->rel_tol_gap && curr_node_idx >= 0 && prob->iter < (prob->stgs->maxit - 1);
}

static int get_ret_code(ecos_bb_pwork *prob)
{
    if (prob->iter < prob->stgs->maxit - 1)
    {
        if (isinf(prob->global_U))
        {
            if (prob->global_U >= 0)
            {
                return MI_INFEASIBLE;
            }
            else
            {
                return MI_UNBOUNDED;
            }
        }
        else
            return MI_OPTIMAL_SOLN;
    }
    else
    {
        if (isinf(prob->global_U))
        {
            if (prob->global_U >= 0)
            {
                return MI_MAXITER_NO_SOLN;
            }
            else
            {
                return MI_MAXITER_UNBOUNDED;
            }
        }
        else
            return MI_MAXITER_FEASIBLE_SOLN;
    }
}

static void initialize_root(ecos_bb_pwork *prob)
{
    idxint i;
    prob->nodes[0].status = MI_NOT_SOLVED;
    prob->nodes[0].L = -ECOS_INFINITY;
    prob->nodes[0].U = ECOS_INFINITY;
    prob->nodes[0].prev_split_idx = -1;
    prob->global_L = -ECOS_INFINITY;
    prob->global_U = ECOS_INFINITY;
    for (i = 0; i < prob->num_bool_vars; ++i)
        prob->bool_node_ids[i] = MI_STAR;
    for (i = 0; i < prob->num_int_vars; ++i)
    {
        prob->int_node_ids[2 * i] = MAX_FLOAT_INT;
        prob->int_node_ids[2 * i + 1] = MAX_FLOAT_INT;
    }
}

idxint ECOS_BB_solve(ecos_bb_pwork *prob)
{
    idxint curr_node_idx = 0;

#if MI_PRINTLEVEL > 0
    if (prob->stgs->verbose)
    {
        PRINTTEXT("Iter\tLower Bound\tUpper Bound\tGap\n");
        PRINTTEXT("================================================\n");
    }
#endif

    /* Initialize to root node and execute steps 1 on slide 6 */
    /* of http://stanford.edu/class/ee364b/lectures/bb_slides.pdf*/
    prob->iter = 0;
    prob->dive_node_id = 0;
    initialize_root(prob);
    /*print_node(prob, curr_node_idx);*/
    get_bounds(curr_node_idx, prob);

    prob->global_L = prob->nodes[curr_node_idx].L;
    prob->global_U = prob->nodes[curr_node_idx].U;

    while (should_continue(prob, curr_node_idx))
    {

#if MI_PRINTLEVEL > 0
        if (prob->stgs->verbose)
        {
            print_progress(prob);
        }
#endif

        ++(prob->iter);

        /* Step 2*/
        /* Branch replaces nodes[curr_node_idx] with leftNode*/
        /* and nodes[prob->iter] with rightNode */
        branch(curr_node_idx, prob);

        /* Step 3*/
        get_bounds(curr_node_idx, prob);
        get_bounds(prob->iter, prob);

        /* Step 4*/
        prob->global_L = get_global_L(prob);

        curr_node_idx = get_next_node(prob);
    }
    load_solution(prob);

#if MI_PRINTLEVEL > 0
    if (prob->stgs->verbose)
    {
        print_progress(prob);
    }
#endif
    return get_ret_code(prob);
}

void updateDataEntry_h(ecos_bb_pwork *w, idxint idx, pfloat value)
{
    ecos_updateDataEntry_h(w->ecos_prob, idx + (2 * w->num_bool_vars), value);
}

void updateDataEntry_c(ecos_bb_pwork *w, idxint idx, pfloat value)
{
    ecos_updateDataEntry_c(w->ecos_prob, idx, value);
}
