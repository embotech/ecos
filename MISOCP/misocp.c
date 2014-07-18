#include "ecos.h"
#include "misocp.h"

// Boolean vars only

void misocp_solve(misocp_pwork* prob){
    prob->iter = 0;

    node* curr_node;
    pfloat* curr_node_id;

    // Initialize to root node and execute steps 1 on slide 6
    idxint curr_node_idx = 0;
    get_bounds(curr_node_idx, prob);

    prob->global_L = prob->nodes[curr_node_idx].L;
    prob->global_U = prob->nodes[curr_node_idx].U;
    
    while ( should_continue(prob) ){
        ++(prob->iter);

        // Step 2
        // Branch replaces nodes[curr_node_idx]=leftNode 
        // and nodes[prob->iter]=rightNode 
        branch(curr_node_idx, prob);

        // Step 3
        get_bounds(curr_node_idx, prob);
        get_bounds(prob->iter, prob);

        // Step 4
        prob->global_L = get_global_L(prob);
        prob->global_U = min(prob->global_U, min(prob->nodes[curr_node_idx].U, prob->nodes[iter].U));

        curr_node_idx = get_next_node(prob);        
    }
}

void branch(idxint curr_node_idx, misocp_pwork* prob){
    
}

/*
 * Function to return the next node to be expanded
 */
idxint get_next_node(misocp_pwork* prob){
    idxint i, next_node;
    pfloat L = INFINITY;
    for(i=0; i < prob->iter; ++i){
        if(nodes[i].status == MI_SOLVED_BRANCHABLE && nodes[i].L < L ){
            next_node = i;
            L = nodes[i].L;
        }
    }
    return next_node;
}

pfloat get_global_L(misocp_pwork* prob){
    idxint i;
    pfloat L = INFINITY;
    for(i=0; i < prob->iter; ++i){ L = min(L, prob->nodes[i].L); }
    return L;
}

/*
 * Function to return the next var to split on
 */
idxint get_branch_var(pfloat* x, idxint num_int_vars){
    idxint i, split_var;
    pfloat d = 1.0;
    for (i=0; i<num_int_vars; ++i){
        pfloat ambiguity = abs(x[i]-0.5);
        if ( ambiguity < d){
            split_var = i;
            d = ambiguity;
        }
    }
    return split_var;
}

void get_bounds(idxint node_idx, misocp_pwork* prob){  
    idxint i, ret_code, branchable;
    set_prob(prob->ecos_prob, get_node_id(node_idx), prob->num_int_vars);
    ret_code = ECOS_solve(ecos_prob);

    branchable = 1;
    if (ret_code == ECOS_OPTIMAL){
        prob->nodes[node_idx].L = prob->ecos_prob->cx;

        // Figure out if x is already an integer solution
        for (i=0; i<prob->num_int_vars; ++i){
            prob->tmp_node_id[i] = (char) round(get_node_id(node_idx)[i]);
            branchable = branchable && 
                !float_eqls(get_node_id(node_idx)[i], (pfloat) prob->tmp_node_id[i]);
        }
        
        if (branchable){ // Round and check feasibility
            prob->nodes[node_idx].split_idx = get_branch_var(get_node_id(node_idx, prob), prob);
            prob->nodes[node_idx].status = MI_SOLVED_BRANCHABLE;  
            set_prob(prob->ecos_prob, prob->tmp_node_id, prob->num_int_vars);
            ret_code = ECOS_solve(ecos_prob);

            if (ret_code == ECOS_OPTIMAL){
                prob->nodes[node_idx].U = prob->ecos_prob->cx;    
            }else{
                prob->nodes[node_idx].U = INFINITY;    
            }
        }else{ // This is already an integer solution
            prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
            prob->nodes[node_idx].U = prob->ecos_prob->cx;
        }

        if (prob->nodes[node_idx].U < prob->global_U){
            for(i=0; i<prob->ecos_prob->n; ++i){ prob->best_x[i] = get_node_id(node_idx)[i]; }
        }
    }else { //Assume node infeasible
        prob->nodes[node_idx].L = INFINITY;
        prob->nodes[node_idx].U = INFINITY;
        prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
    }
}

idxint should_continue(misocp_pwork* prob){
    return (prob->global_U - prob->global_L) > MI_ABS_EPS 
        || (prob->global_U / prob->global_L - 1.0) > MI_REL_EPS
        || prob->iter < MI_MAXITER;
}

/* 
 * Updates the solver's lb and ub contraints for integer variables
 * to the bounds specified by the node
 */
void set_prob(pwork* ecos_prob, char* node_id, idxint num_int_vars){
    for(i=0; i<num_int_vars; ++i){
        switch(node_id[i]){
            case ONE:
                ecos_prob->b[2*i] = -((pfloat) ONE); // -x <= -1 <=> x >= 1
                ecos_prob->b[2*i + 1] = (pfloat) ONE;         
                break;
            case ZERO:
                ecos_prob->b[2*i] = (pfloat) ZERO;
                ecos_prob->b[2*i + 1] = (pfloat) ZERO;         
                break;
            case STAR:
                ecos_prob->b[2*i] = (pfloat) ZERO;
                ecos_prob->b[2*i + 1] = (pfloat) ONE;         
                break;
        }
    }
}
