#include "glblopts.h"
#include "ecos.h"
#include "ecos_bb.h"
#include "math.h"
#include "equil.h"
#include "spla.h"

// Boolean vars only

void branch(idxint curr_node_idx, ecos_bb_pwork* prob){
    idxint i;

    // Create right node
    prob->nodes[prob->iter].L = prob->nodes[curr_node_idx].L;
    prob->nodes[prob->iter].status = MI_NOT_SOLVED;

    // Copy over the node id
    for(i=0; i < prob->num_bool_vars; ++i){
        get_node_id(prob->iter, prob)[i] = get_node_id(curr_node_idx, prob)[i];
    }
    
    get_node_id(curr_node_idx, prob)[prob->nodes[curr_node_idx].split_idx] = MI_ZERO;
    get_node_id(prob->iter, prob)[prob->nodes[curr_node_idx].split_idx] = MI_ONE;
    
    prob->nodes[curr_node_idx].status = MI_NOT_SOLVED;
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
    for(i=0; i <= prob->iter; ++i) L = min(L,prob->nodes[i].L);
    return L;
}

/*
 * Function to return the next var to split on
 */
idxint get_branch_var(pfloat* x, idxint num_bool_vars){
    idxint i, split_var;
    pfloat d = 1.0;
    for (i=0; i<num_bool_vars; ++i){
        pfloat ambiguity = abs_2(x[i]-0.5);
        if ( ambiguity < d){
            split_var = i;
            d = ambiguity;
        }
    }
    return split_var;
}

/* 
 * Updates the solver's lb and ub contraints for integer variables
 * to the bounds specified by the node
 */
void set_prob(pwork* ecos_prob, char* node_id, idxint num_bool_vars){
    idxint i;
    //unset_equilibration(ecos_prob);
    for(i=0; i<num_bool_vars; ++i){
        switch(node_id[i]){
            case MI_ONE:
                ecos_updateDataEntry_h(ecos_prob, 2*i, -1.0); //ecos_prob->h[2*i] = -(1.0); // -x <= -1 <=> x >= 1
                ecos_updateDataEntry_h(ecos_prob, 2*i+1, 1.0);//ecos_prob->h[2*i + 1] = 1.0;
                break;
            case MI_ZERO:
                ecos_updateDataEntry_h(ecos_prob, 2*i, 0.0);//ecos_prob->h[2*i] = 0.0;
                ecos_updateDataEntry_h(ecos_prob, 2*i+1, 0.0);//ecos_prob->h[2*i + 1] = 0.0;
                break;
            case MI_STAR:
                ecos_updateDataEntry_h(ecos_prob, 2*i, 0.0);//ecos_prob->h[2*i] = 0.0;
                ecos_updateDataEntry_h(ecos_prob, 2*i+1, 1.0);//ecos_prob->h[2*i + 1] = 1.0;
                break;
        }
        //PRINTTEXT("lb:%f, ub:%f\n", ecos_prob->h[2*i], ecos_prob->h[2*i+1]);
    }
    //set_equilibration(ecos_prob);
}

void get_bounds(idxint node_idx, ecos_bb_pwork* prob){  
    idxint i, ret_code, branchable;
    set_prob(prob->ecos_prob, get_node_id(node_idx,prob), prob->num_bool_vars);    
    ret_code = ECOS_solve(prob->ecos_prob);

    
    if (ret_code == ECOS_OPTIMAL){
        prob->nodes[node_idx].L = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);

        branchable = 1;
        // Figure out if x is already an integer solution
        for (i=0; i<prob->num_bool_vars; ++i){
            prob->tmp_node_id[i] = (char) round( prob->ecos_prob->x[i] );
            branchable &= float_eqls( prob->ecos_prob->x[i] , (pfloat) prob->tmp_node_id[i]);
            //PRINTTEXT("branch: %f %f\n", prob->ecos_prob->x[i] , (pfloat) prob->tmp_node_id[i]);
        }
        branchable = !branchable;
        //PRINTTEXT("branchable: %u\n", branchable);

        //PRINTTEXT("Orig Solve: %u\n", ret_code);for (i=0; i<prob->ecos_prob->n; ++i) PRINTTEXT("%f\n", prob->ecos_prob->x[i]);

        if (branchable){ // Round and check feasibility
            prob->nodes[node_idx].split_idx = get_branch_var(prob->ecos_prob->x, prob->num_bool_vars);
            prob->nodes[node_idx].status = MI_SOLVED_BRANCHABLE;  
            set_prob(prob->ecos_prob, prob->tmp_node_id, prob->num_bool_vars);
            ret_code = ECOS_solve(prob->ecos_prob);

            //PRINTTEXT("Guess: %u\n", ret_code);for (i=0; i<prob->ecos_prob->n; ++i) PRINTTEXT("%f\n", prob->ecos_prob->x[i]);

            if (ret_code == ECOS_OPTIMAL){
                prob->nodes[node_idx].U = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);
            }else{
                prob->nodes[node_idx].U = INFINITY;    
            }
        }else{ // This is already an integer solution
            prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
            prob->nodes[node_idx].U = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);
        }

        if (prob->nodes[node_idx].U < prob->global_U){
            for(i=0; i<prob->ecos_prob->n; ++i){ prob->best_x[i] = prob->ecos_prob->x[i]; }
        }

    }else { //Assume node infeasible
        //PRINTTEXT("Ret code: %u\n", ret_code);

        prob->nodes[node_idx].L = INFINITY;
        prob->nodes[node_idx].U = INFINITY;
        prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
    }

    //PRINTTEXT("%f, %f, %f, %f\n", prob->nodes[node_idx].L, prob->nodes[node_idx].U, prob->nodes[node_idx].U - prob->nodes[node_idx].L, prob->nodes[node_idx].U/prob->nodes[node_idx].L-1.0);
}

idxint should_continue(ecos_bb_pwork* prob, idxint curr_node_idx){
    return (prob->global_U - prob->global_L) > MI_ABS_EPS 
        && abs_2(prob->global_U / prob->global_L - 1.0) > MI_REL_EPS
        && curr_node_idx >= 0
        && prob->iter < MI_MAXITER;
}

void print_progress(ecos_bb_pwork* prob){
    PRINTTEXT("%u \t%.2f \t\t%.2f \t\t%.2f\n", (int)prob->iter, prob->global_L, prob->global_U, prob->global_U-prob->global_L);
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


void print_node(ecos_bb_pwork* prob, idxint i){
    int j;
    PRINTTEXT("Node info: %u : %f : %f : %u\nPartial id:", prob->nodes[i].status, 
        prob->nodes[i].L, prob->nodes[i].U, (int)prob->nodes[i].split_idx);
    for (j=0; j<prob->num_bool_vars; ++j) PRINTTEXT("%i,", get_node_id(i,prob)[j]);
    PRINTTEXT("\n");
}

void initialize_root(ecos_bb_pwork* prob){
    idxint i;
    prob->nodes[0].status = MI_NOT_SOLVED;
    prob->nodes[0].L = -INFINITY;
    prob->nodes[0].U =  INFINITY;
    for (i=0; i < prob->num_bool_vars; ++i){ prob->node_ids[i] = MI_STAR; }
}

int ecos_bb_solve(ecos_bb_pwork* prob){
    
    prob->iter = 0;
    
    // Initialize to root node and execute steps 1 on slide 6 
    // of http://stanford.edu/class/ee364b/lectures/bb_slides.pdf
    initialize_root(prob);
    idxint curr_node_idx = 0;
    //print_node(prob, curr_node_idx);
    get_bounds(curr_node_idx, prob);

    prob->global_L = prob->nodes[curr_node_idx].L;
    prob->global_U = prob->nodes[curr_node_idx].U;

    PRINTTEXT("Iter\tLower Bound\tUpper Bound\tGap\n");
    PRINTTEXT("================================================\n");
    while ( should_continue(prob, curr_node_idx) ){
        print_progress(prob);

        ++(prob->iter);

        // Step 2
        // Branch replaces nodes[curr_node_idx] with leftNode
        // and nodes[prob->iter] with rightNode 
        branch(curr_node_idx, prob);

        //PRINTTEXT("curr_node_idx: %u\n", curr_node_idx);

        // Step 3
        get_bounds(curr_node_idx, prob);
        get_bounds(prob->iter, prob);

        //for (i=0; i <= prob->iter; ++i) print_node(prob, i);

        // Step 4
        prob->global_L = get_global_L(prob);
        prob->global_U = min(prob->global_U, min(prob->nodes[curr_node_idx].U, prob->nodes[prob->iter].U));

        curr_node_idx = get_next_node(prob);        
    }
    print_progress(prob);
    return get_ret_code(prob);
}


