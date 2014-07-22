#include "ecos.h"
#include "misocp.h"
#include "math.h"
#include "equil.h"
#include "spla.h"

// Boolean vars only

void branch(idxint curr_node_idx, misocp_pwork* prob){
    idxint i;

    // Create right node
    prob->nodes[prob->iter].L = prob->nodes[curr_node_idx].L;
    prob->nodes[prob->iter].status = MI_NOT_SOLVED;

    // Copy over the node id
    for(i=0; i<prob->num_int_vars; ++i){
        get_node_id(prob->iter, prob)[i] = get_node_id(curr_node_idx, prob)[i];
    }

    get_node_id(curr_node_idx, prob)[prob->nodes[curr_node_idx].split_idx] = MI_ZERO;
    get_node_id(prob->iter, prob)[prob->nodes[curr_node_idx].split_idx] = MI_ONE;

    prob->nodes[curr_node_idx].status = MI_NOT_SOLVED;
}   

/*
 * Function to return the next node to be expanded
 */
idxint get_next_node(misocp_pwork* prob){
    idxint i;
    idxint next_node = -1;
    pfloat L = INFINITY;
    for(i=0; i < prob->iter; ++i){
        if(prob->nodes[i].status == MI_SOLVED_BRANCHABLE && prob->nodes[i].L < L ){
            next_node = i;
            L = prob->nodes[i].L;
        }
    }
    return next_node;
}

pfloat get_global_L(misocp_pwork* prob){
    idxint i;
    pfloat L = INFINITY;
    for(i=0; i <= prob->iter; ++i) L = min(L,prob->nodes[i].L);
    return L;
}

/*
 * Function to return the next var to split on
 */
idxint get_branch_var(pfloat* x, idxint num_int_vars){
    idxint i, split_var;
    pfloat d = 1.0;
    for (i=0; i<num_int_vars; ++i){
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
void set_prob(pwork* ecos_prob, char* node_id, idxint num_int_vars){
    idxint i;
    unset_equilibration(ecos_prob);
    for(i=0; i<num_int_vars; ++i){
        switch(node_id[i]){
            case MI_ONE:
                ecos_prob->h[2*i] = -(1.0); // -x <= -1 <=> x >= 1
                ecos_prob->h[2*i + 1] = 1.0;         
                break;
            case MI_ZERO:
                ecos_prob->h[2*i] = 0.0;
                ecos_prob->h[2*i + 1] = 0.0;         
                break;
            case MI_STAR:
                ecos_prob->h[2*i] = 0.0;
                ecos_prob->h[2*i + 1] = 1.0;         
                break;
        }
    }
    set_equilibration(ecos_prob);
}

void get_bounds(idxint node_idx, misocp_pwork* prob){  
    idxint i, ret_code, branchable;
    set_prob(prob->ecos_prob, get_node_id(node_idx,prob), prob->num_int_vars);
    
    ret_code = ECOS_solve(prob->ecos_prob);



    branchable = 1;
    if (ret_code == ECOS_OPTIMAL){
        prob->nodes[node_idx].L = eddot(prob->ecos_prob->n, prob->ecos_prob->x, prob->ecos_prob->c);

        // Figure out if x is already an integer solution
        for (i=0; i<prob->num_int_vars; ++i){
            prob->tmp_node_id[i] = (char) round( prob->ecos_prob->x[i] );
            branchable = branchable && 
                !float_eqls( prob->ecos_prob->x[i] , (pfloat) prob->tmp_node_id[i]);
        }
        
       // printf("\n");for (i=0; i<prob->ecos_prob->n; ++i) printf("%f\n", prob->ecos_prob->x[i]);

        if (branchable){ // Round and check feasibility
            prob->nodes[node_idx].split_idx = get_branch_var(prob->ecos_prob->best_x, prob->num_int_vars);
            prob->nodes[node_idx].status = MI_SOLVED_BRANCHABLE;  
            set_prob(prob->ecos_prob, prob->tmp_node_id, prob->num_int_vars);
            ret_code = ECOS_solve(prob->ecos_prob);

            //printf("\n");for (i=0; i<prob->ecos_prob->n; ++i) printf("%f\n", prob->ecos_prob->x[i]);

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
        prob->nodes[node_idx].L = INFINITY;
        prob->nodes[node_idx].U = INFINITY;
        prob->nodes[node_idx].status = MI_SOLVED_NON_BRANCHABLE;
    }

    //printf("%f, %f, %f, %f\n", prob->nodes[node_idx].L, prob->nodes[node_idx].U, prob->nodes[node_idx].U - prob->nodes[node_idx].L, prob->nodes[node_idx].U/prob->nodes[node_idx].L-1.0);
}

idxint should_continue(misocp_pwork* prob, idxint curr_node_idx){
    return (prob->global_U - prob->global_L) > MI_ABS_EPS 
        && abs_2(prob->global_U / prob->global_L - 1.0) > MI_REL_EPS
        && curr_node_idx >= 0
        && prob->iter < MI_MAXITER;
}

void print_progress(misocp_pwork* prob){
    printf("iter: %u, L: %f, U: %f\n", prob->iter, prob->global_L, prob->global_U);
}

int get_ret_code(misocp_pwork* prob){
    if ( isinf(prob->global_U) ){ return ECOS_OPTIMAL;
    }else if ( prob->iter == MI_MAXITER ){ return ECOS_MAXIT;        
    }else { return ECOS_OPTIMAL; }
}

void print_node(misocp_pwork* prob, idxint i){
    int j;
    printf("\n%u : %f : %f : %u\n", prob->nodes[i].status, 
        prob->nodes[i].L, prob->nodes[i].U, prob->nodes[i].split_idx);
    for (j=0; j<prob->num_int_vars; ++j) printf("%i,", get_node_id(i,prob)[j]);
    printf("\n\n");
}

int misocp_solve(misocp_pwork* prob){
    prob->iter = 0;
    
    // Initialize to root node and execute steps 1 on slide 6
    idxint curr_node_idx = 0;
    get_bounds(curr_node_idx, prob);

    prob->global_L = prob->nodes[curr_node_idx].L;
    prob->global_U = prob->nodes[curr_node_idx].U;
    
    print_node(prob, 0);

    while ( should_continue(prob, curr_node_idx) ){
        print_progress(prob);

        ++(prob->iter);

        // Step 2
        // Branch replaces nodes[curr_node_idx]=leftNode 
        // and nodes[prob->iter]=rightNode 
        branch(curr_node_idx, prob);

        print_node(prob, 0);
        print_node(prob,1);

        // Step 3
        get_bounds(curr_node_idx, prob);
        get_bounds(prob->iter, prob);

        print_node(prob, 0);
        print_node(prob,1);

        // Step 4
        prob->global_L = get_global_L(prob);
        prob->global_U = min(prob->global_U, min(prob->nodes[curr_node_idx].U, prob->nodes[prob->iter].U));

        curr_node_idx = get_next_node(prob);        
    }
    print_progress(prob);
    return get_ret_code(prob);
}


