#include "glblopts.h"
#include "ecos.h"
#include "ecos_bb.h"
#include "math.h"
#include "equil.h"
#include "spla.h"

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

/* Print utility functions*/
#if MI_PRINTLEVEL > 0
void print_progress(ecos_bb_pwork* prob){
    PRINTTEXT("%u \t%.2f \t\t%.2f \t\t%.2f\n", (int)prob->iter, prob->global_L, prob->global_U, prob->global_U-prob->global_L);
}

void print_ecos_solution(ecos_bb_pwork* prob){
    int i; PRINTTEXT("ecos->x: ");
    for(i=0; i < prob->ecos_prob->n; ++i) PRINTTEXT("%.2f ", prob->ecos_prob->x[i] );
    PRINTTEXT("\n");
}

void print_ecos_xequil(ecos_bb_pwork* prob){
    int i; PRINTTEXT("ecos->xequil: ");
    for (i=0; i<prob->ecos_prob->n; ++i) PRINTTEXT("%.2f ", prob->ecos_prob->xequil[i] );
    PRINTTEXT("\n");
}


void print_ecos_h(ecos_bb_pwork* prob){
    int i; PRINTTEXT("ecos->h: ");
    for (i=0; i<prob->ecos_prob->m; ++i) PRINTTEXT("%.2f ", prob->ecos_prob->h[i] );
    PRINTTEXT("\n");
}

void print_ecos_c(ecos_bb_pwork* prob){
    int i; PRINTTEXT("ecos->c: ");
    for (i=0; i<prob->ecos_prob->n; ++i) PRINTTEXT("%.2f ", prob->ecos_prob->c[i] );
    PRINTTEXT("\n");
}

void print_node(ecos_bb_pwork* prob, idxint i){
    int j; PRINTTEXT("Node info: %u : %.2f : %.2f : %u : %.2f Partial id:", prob->nodes[i].status,
        prob->nodes[i].L, prob->nodes[i].U, (int)prob->nodes[i].split_idx, prob->nodes[i].split_val);
    for(j=0; j < prob->num_bool_vars; ++j)
        PRINTTEXT("%i ", get_bool_node_id(i, prob)[j] );
    PRINTTEXT(" | ");
    for(j=0; j < prob->num_int_vars; ++j)
        PRINTTEXT("(%.2f, %.2f) ", -get_int_node_id(i, prob)[2*j], get_int_node_id(i, prob)[2*j+1] );
    PRINTTEXT("\n");
}
#endif


void branch(idxint curr_node_idx, ecos_bb_pwork* prob){
    idxint i, split_idx = prob->nodes[curr_node_idx].split_idx;

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose) {
        PRINTTEXT("Branching->\t");
        print_node(prob, curr_node_idx);
    }
#endif

    /* Create right node*/
    prob->nodes[prob->iter].L = prob->nodes[curr_node_idx].L;
    prob->nodes[prob->iter].status = MI_NOT_SOLVED;

    /* Copy over the node id*/
    for(i=0; i < prob->num_bool_vars; ++i)
        get_bool_node_id(prob->iter, prob)[i] = get_bool_node_id(curr_node_idx, prob)[i];
    for(i=0; i < prob->num_int_vars*2; ++i)
        get_int_node_id(prob->iter, prob)[i] = get_int_node_id(curr_node_idx, prob)[i];

    if (split_idx < prob->num_bool_vars){
        get_bool_node_id(curr_node_idx, prob)[split_idx] = MI_ZERO;
        get_bool_node_id(prob->iter, prob)[split_idx] = MI_ONE;
    }else{
        split_idx -= prob->num_bool_vars;

        /* Left branch constrain UB */
        get_int_node_id(curr_node_idx, prob)[split_idx*2 + 1] =
            pfloat_floor( prob->nodes[curr_node_idx].split_val );

        /* Right branch constrain LB */
        get_int_node_id(prob->iter, prob)[split_idx*2 ] =
            -pfloat_ceil( prob->nodes[curr_node_idx].split_val );
    }

    prob->nodes[curr_node_idx].status = MI_NOT_SOLVED;

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose) {
        PRINTTEXT(" Left-> \t "); print_node(prob, curr_node_idx);
        PRINTTEXT(" Right->\t "); print_node(prob, prob->iter);
    }
#endif
}

/*
 * Function to return the next node to be expanded
 */
idxint get_next_node(ecos_bb_pwork* prob){
    idxint i;
    idxint next_node = -1;
    pfloat L = INFINITY;
    for(i=0; i <= prob->iter; ++i){
        if(prob->nodes[i].status == MI_SOLVED_BRANCHABLE && prob->nodes[i].L < L ){
            next_node = i;
            L = prob->nodes[i].L;
        }
    }
    return next_node;
}

pfloat get_global_L(ecos_bb_pwork* prob){
    idxint i;
    pfloat L = INFINITY;
    for(i=0; i <= prob->iter; ++i) L = MIN(L,prob->nodes[i].L);
    return L;
}

/*
 * Function to return the next var to split on
 */
void get_branch_var(ecos_bb_pwork* prob, idxint* split_idx, pfloat* split_val){
    idxint i;
    pfloat x, y, d, ambiguity;
    d = 1.0;
    for (i=0; i<(prob->num_bool_vars + prob->num_int_vars); ++i){
        if (i < prob->num_bool_vars ){
            y = prob->ecos_prob->x[ prob->bool_vars_idx[i] ];
            x = y;
        } else {
            y = prob->ecos_prob->x[ prob->int_vars_idx[i] ];
            x = y - pfloat_floor(y);
        }
        ambiguity = abs_2(x-0.5);
        if ( ambiguity < d){
            *split_idx = i;
            *split_val = y;
            d = ambiguity;
        }
    }
}

/*
 * Updates the solver's lb and ub contraints for integer variables
 * to the bounds specified by the node
 */
void set_prob(ecos_bb_pwork* prob, char* bool_node_id, pfloat* int_node_id){
    idxint i;
    for(i=0; i<prob->num_bool_vars; ++i){
        switch(bool_node_id[i]){
            case MI_ONE:
                ecos_updateDataEntry_h(prob->ecos_prob, 2*i, -1.0); /*ecos_prob->h[2*i] = -(1.0); // -x <= -1 <=> x >= 1*/
                ecos_updateDataEntry_h(prob->ecos_prob, 2*i+1, 1.0);/*ecos_prob->h[2*i + 1] = 1.0;*/
                break;
            case MI_ZERO:
                ecos_updateDataEntry_h(prob->ecos_prob, 2*i, 0.0);/*ecos_prob->h[2*i] = 0.0;*/
                ecos_updateDataEntry_h(prob->ecos_prob, 2*i+1, 0.0);/*ecos_prob->h[2*i + 1] = 0.0;*/
                break;
            case MI_STAR:
                ecos_updateDataEntry_h(prob->ecos_prob, 2*i, 0.0);/*ecos_prob->h[2*i] = 0.0;*/
                ecos_updateDataEntry_h(prob->ecos_prob, 2*i+1, 1.0);/*ecos_prob->h[2*i + 1] = 1.0;*/
                break;
			default:
#if MI_PRINTLEVEL > 1
				PRINTTEXT("Illegal boolean setting arguments passed: %u \n", bool_node_id[i]);
#endif
				break;
        }
    }

    /* Set integer bounds */
    for(i=0; i<prob->num_int_vars; ++i){
        ecos_updateDataEntry_h(prob->ecos_prob, 2*(i + prob->num_bool_vars) , int_node_id[2*i] );
        ecos_updateDataEntry_h(prob->ecos_prob, 2*(i + prob->num_bool_vars)+1 , int_node_id[2*i+1] );
    }

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose){ print_ecos_h(prob); }
#endif

}

/*
 * Stores the ecos solution to the array inside ecos_bb
 */
void store_solution(ecos_bb_pwork* prob){
    idxint i;
    for(i=0; i<prob->ecos_prob->n; ++i) prob->x[i] = prob->ecos_prob->x[i];
    for(i=0; i<prob->ecos_prob->p; ++i) prob->y[i] = prob->ecos_prob->y[i];
    for(i=0; i<prob->ecos_prob->m; ++i) prob->z[i] = prob->ecos_prob->z[i];
    for(i=0; i<prob->ecos_prob->m; ++i) prob->s[i] = prob->ecos_prob->s[i];
    prob->kap  =  prob->ecos_prob->kap ;
    prob->tau  =  prob->ecos_prob->tau ;
    *(prob->best_info) = *(prob->ecos_prob->best_info);
}

/*
 * Loads the ecos_bb solution back into the ecos struct, necessary for the Matlab/Python interface
 */
void load_solution(ecos_bb_pwork* prob){
    idxint i;
    for(i=0; i<prob->ecos_prob->n; ++i) prob->ecos_prob->x[i] = prob->x[i];
    for(i=0; i<prob->ecos_prob->p; ++i) prob->ecos_prob->y[i] = prob->y[i];
    for(i=0; i<prob->ecos_prob->m; ++i) prob->ecos_prob->z[i] = prob->z[i];
    for(i=0; i<prob->ecos_prob->m; ++i) prob->ecos_prob->s[i] = prob->s[i];
    prob->ecos_prob->kap = prob->kap;
    prob->ecos_prob->tau = prob->tau;
    *(prob->ecos_prob->best_info) = *(prob->best_info);
}

void get_bounds(idxint node_idx, ecos_bb_pwork* prob){
    idxint i, ret_code, branchable;
    set_prob( prob, get_bool_node_id(node_idx,prob), get_int_node_id(node_idx, prob) );
    ret_code = ECOS_solve(prob->ecos_prob);

#if MI_PRINTLEVEL > 1
    if (prob->stgs->verbose){ print_ecos_solution(prob); }
    if (ret_code != ECOS_OPTIMAL && ret_code != ECOS_PINF && ret_code != ECOS_DINF){
        PRINTTEXT("Exit code: %u\n", ret_code);
    }
#endif


    if (ret_code == ECOS_OPTIMAL || 
        ret_code == ECOS_INACC_OFFSET || 
        ret_code == ECOS_MAXIT ||
        ret_code == ECOS_NUMERICS ){
        prob->nodes[node_idx].L = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);

        /* Figure out if x is already an integer solution if the solution had no numerical errors*/
        branchable = 1;
        if (ret_code == ECOS_OPTIMAL){
            for (i=0; i<prob->num_bool_vars; ++i){
                prob->tmp_bool_node_id[i] = (char) pfloat_round( prob->ecos_prob->x[prob->bool_vars_idx[i]] );
                branchable &= float_eqls( prob->ecos_prob->x[i] , (pfloat) prob->tmp_bool_node_id[i] );
            }
            for (i=0; i<prob->num_int_vars; ++i){
				prob->tmp_int_node_id[2 * i + 1] = pfloat_round(prob->ecos_prob->x[prob->int_vars_idx[i]]);
                prob->tmp_int_node_id[2*i] = -(prob->tmp_int_node_id[2*i + 1]);
                branchable &= float_eqls( prob->ecos_prob->x[i] , prob->tmp_int_node_id[2*i + 1] );
            }
            branchable = !branchable;
        }

        if (branchable){ /* pfloat_round and check feasibility*/
            get_branch_var(prob, &(prob->nodes[node_idx].split_idx), &(prob->nodes[node_idx].split_val) );
            prob->nodes[node_idx].status = MI_SOLVED_BRANCHABLE;

#if MI_PRINTLEVEL > 1
            if (prob->stgs->verbose){ PRINTTEXT("Rounded Solution:\n"); }
#endif

            set_prob(prob, prob->tmp_bool_node_id, prob->tmp_int_node_id);
            ret_code = ECOS_solve(prob->ecos_prob);

#if MI_PRINTLEVEL > 1
            if (prob->stgs->verbose){ print_ecos_solution(prob); }
            if (ret_code != ECOS_OPTIMAL && ret_code != ECOS_PINF && ret_code != ECOS_OPTIMAL){
                PRINTTEXT("Exit code: %u\n", ret_code);
            }
#endif

            if (ret_code == ECOS_OPTIMAL){
                /* Use the node's U as tmp storage */
                prob->nodes[node_idx].U = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);

#if MI_PRINTLEVEL > 1
                if (prob->stgs->verbose){
                    PRINTTEXT("New optimal solution, U: %.2f\n", prob->nodes[node_idx].U);
                    print_ecos_xequil(prob);
                    print_ecos_c(prob);
                    print_ecos_solution(prob);
                }
#endif

                store_solution(prob);
                prob->global_U = prob->nodes[node_idx].U;
            }

            /* Reset the node's U back to INF because it was not originally feasible */
            prob->nodes[node_idx].U = INFINITY;
        }else{ /* This is already an integer solution*/
            prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
            prob->nodes[node_idx].U = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);
        }

        if (prob->nodes[node_idx].U < prob->global_U){

#if MI_PRINTLEVEL > 1
            if (prob->stgs->verbose){
                PRINTTEXT("New optimal solution, U: %.2f\n", prob->nodes[node_idx].U);
                print_ecos_xequil(prob);
                print_ecos_c(prob);
                print_ecos_solution(prob);
            }
#endif

            store_solution(prob);
            prob->global_U = prob->nodes[node_idx].U;
        }

    }else { /*Assume node infeasible*/
        prob->nodes[node_idx].L = INFINITY;
        prob->nodes[node_idx].U = INFINITY;
        prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
    }
}

idxint should_continue(ecos_bb_pwork* prob, idxint curr_node_idx){
    return (prob->global_U - prob->global_L) > MI_ABS_EPS
        && abs_2(prob->global_U / prob->global_L - 1.0) > MI_REL_EPS
        && curr_node_idx >= 0
        && prob->iter < (prob->maxiter-1);
}

int get_ret_code(ecos_bb_pwork* prob){
    if ( prob->iter < MI_MAXITER){
        if ( isinf(prob->global_U) ) return MI_INFEASIBLE;
        else return MI_OPTIMAL_SOLN;
    } else {
        if ( isinf(prob->global_U) ) return MI_MAXITER_NO_SOLN;
        else return MI_MAXITER_FEASIBLE_SOLN;
    }
}

void initialize_root(ecos_bb_pwork* prob){
    idxint i;
    prob->nodes[0].status = MI_NOT_SOLVED;
    prob->nodes[0].L = -INFINITY;
    prob->nodes[0].U =  INFINITY;
    prob->global_L = -INFINITY;
    prob->global_U = INFINITY;
    for (i=0; i < prob->num_bool_vars; ++i) prob->bool_node_ids[i] = MI_STAR;
    for (i=0; i < prob->num_int_vars; ++i) {prob->int_node_ids[2*i] = MAX_FLOAT_INT; prob->int_node_ids[2*i+1] = MAX_FLOAT_INT;}
}

idxint ECOS_BB_solve(ecos_bb_pwork* prob) {
    idxint curr_node_idx = 0;

#if MI_PRINTLEVEL > 0
    if (prob->stgs->verbose){
        PRINTTEXT("Iter\tLower Bound\tUpper Bound\tGap\n");
        PRINTTEXT("================================================\n");
    }
#endif

    /* Initialize to root node and execute steps 1 on slide 6 */
    /* of http://stanford.edu/class/ee364b/lectures/bb_slides.pdf*/
    prob->iter = 0;
    initialize_root(prob);
    /*print_node(prob, curr_node_idx);*/
    get_bounds(curr_node_idx, prob);

    prob->global_L = prob->nodes[curr_node_idx].L;
    prob->global_U = prob->nodes[curr_node_idx].U;

    while ( should_continue(prob, curr_node_idx) ){

#if MI_PRINTLEVEL > 0
        if (prob->stgs->verbose){ print_progress(prob); }
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
    if (prob->stgs->verbose){ print_progress(prob); }
#endif

    return get_ret_code(prob);
}

void updateDataEntry_h(ecos_bb_pwork* w, idxint idx, pfloat value){
    ecos_updateDataEntry_h(w->ecos_prob, idx + (2*w->num_bool_vars), value);
}

void updateDataEntry_c(ecos_bb_pwork* w, idxint idx, pfloat value){
    ecos_updateDataEntry_c(w->ecos_prob, idx , value);
}
