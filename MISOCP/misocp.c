#include "ecos.h"
#include "stdlib.h"
#include "math.h"
#include "misocp.h"

// We assume that the first num_int column vars are integers

/* 
 * ASSUMES THAT G HAS ALREADY BEEN PROCESSED BY socp_to_misocp()!!
 *
 * Modify the input parameters before creating the template ECOS problem
 */
misocp_pwork* misocp_setup(
    idxint n, idxint m, idxint p, 
    idxint l, idxint ncones, idxint* q,
    pfloat* Gpr, idxint* Gjc, idxint* Gir,
    pfloat* Apr, idxint* Ajc, idxint* Air,
    pfloat* c, pfloat* h, pfloat* b, idxint num_int)
{
    idxint i;
    
    // Malloc the problem's memory
    misocp_pwork* prob = (misocp_pwork*) malloc(sizeof(misocp_pwork));

    // Setup the ecos solver
    prob->ecos_prob = ECOS_setup(
        n, m, p, l, ncones, q,
        Gpr, Gjc, Gir,
        Apr, Ajc, Air,
        c, h, b);

    // Malloc the initial nodes
    prob->nodes = (misocp_node*) malloc(INITIAL_NODES*sizeof(misocp_node));

    // Initialize the nodes to their default values
    for (i=0; i<INITIAL_NODES; ++i){

        prob->nodes[i].ret_code = EMPTY;
        prob->nodes[i].last_split = 0;

        // Allocate enough space in this node to hold all the constraints
        // for identified integer variables
        prob->nodes[i].constraints = 
            (pfloat*) malloc(2*num_int*sizeof(pfloat));
    }

    // Store the length of the nodes array
    prob->num_nodes = INITIAL_NODES; 

    // Initialize the tmp_node
    prob->tmp_node.constraints = (pfloat*) malloc(2*num_int*sizeof(pfloat));

    // Store the number of integer variables in the problem
    prob->num_int = num_int;

    // Now initialize the root node
    prob->nodes[i].ret_code = VALID;
    prob->nodes[i].lb = -1e999;
    for (i=0; i<num_int; ++i){
        prob->nodes[i].constraints[2*i] = 1e999;
        prob->nodes[i].constraints[2*i + 1] = -1e999;
    }

    return prob;
}

// Performs the same function as ecos_cleanup
void misocp_cleanup(misocp_pwork* prob){
    idxint i;

    // Free solver memory
    ecos_cleanup(prob->ecos_prob);

    // Free constraints array of every node
    for (i=0; i < prob->num_nodes; ++i){
        free(prob->nodes[i].constraints);
    }

    // Free the nodes array
    free(prob->nodes);

    // Free the constraints array for the tmp node
    free(prob->tmp_node.constraints);

    // Free the remainder of the memory
    free(prob);
}


// This function will extend the nodes array of a problem
// by EXPAND_NODES amount
void resize_node_array(misocp_pwork* prob){
    idxint i;

    // Resize the array
    realloc(prob->nodes, 
        (prob->num_nodes + EXPAND_NODES)*sizeof(misocp_node));

    // Initialize the new nodes
    for (i=num_nodes; i < (num_nodes+EXPAND_NODES); ++i){
        
        prob->nodes[i].ret_code = EMPTY;
        prob->nodes[i].last_split = 0;

        // Allocate enough space in this node to hold all the constraints
        // for identified integer variables
        prob->nodes[i].constraints = 
            (pfloat*) malloc(2*(prob->num_int)*sizeof(pfloat));
    }

    prob->num_nodes += EXPAND_NODES;
}


/* Augments the G and b arrays to take lb and ub constraints
 * for all of the variables marked integer
 * USES MALLOC
 */ 
void socp_to_misocp(idxint num_int, idxint n, idxint m, idxint* l,
    pfloat* Gpr_in, idxint* Gjc_in, idxint* Gir_in,
    pfloat* Gpr_out, idxint* Gjc_out, idxint* Gir_out
    pfloat* b_in, pfloat* b_out)
{
    idxint i;
    /** We assume that the first num_int column vars are integers **/
    
    // First we need to resize the G mat to accomodate lb and ub constraints
    idxint new_G_size = Gjc[n] + 2 * num_int;
    Gpr_out = (pfloat *) malloc( new_G_size * sizeof(pfloat) );
    Gjc_out = (idxint *) malloc( n * sizeof(idxint) );
    Gir_out = (idxint *) malloc( new_G_size * sizeof(idxint) );

    b_out = (pfloat *) malloc( (n + 2 * num_int ) * sizeof(pfloat) );

    // Fill in column pointers making room for the lb's and ub's
    for (i=0; i<n ; ++i){
        Gjc_out[i] = Gjc[i] + 2 * min(i, num_int); 
    }

    // Fill in the new matrices 
    for (i=0; i<num_int; ++i){
        if (i<num_int){
            // Set coefficients for G matrix
            Gpr_out[ Gjc_out[i] ] = 1;
            Gpr_out[ Gjc_out[i] + 1 ] = -1;

            Gir_out[ Gjc_out[i] ] = 2*i;
            Gir_out[ Gjc_out[i] + 1 ] = 2*i + 1;

            // Copy data
            memcpy(&Gpr[Gjr[i]] , &Gpr_out[Gjr_out[i] + 2]
                (Gjr[i+1] - Gjr[i]) *sizeof(pfloat) );        

            // Copy indexes
            memcpy(&Gir[Gjr[i]] , &Gir_out[Gjr_out[i] + 2]
                (Gjr[i+1] - Gjr[i]) *sizeof(idxint) );        

            // Set upper bound to +Inf
            b_out[ 2*i ] = 1e999;     

            // Set lower bound to -Inf
            b_out[ 2*i + 1] = -1e999;     
        } else {
            // Copy data
            memcpy(&Gpr[Gjr[i]] , &Gpr_out[ Gjr_out[i] ]
                (Gjr[i+1] - Gjr[i]) *sizeof(pfloat) );        

            // Copy indexes
            memcpy(&Gir[Gjr[i]] , &Gir_out[ Gjr_out[i] ]
                (Gjr[i+1] - Gjr[i]) *sizeof(idxint) );        

        }
    }

    // Copy the remaining entries of b
    memcpy( b_in, &b_out[2*num_int], m*sizeof(pfloat) );

    // Augment the number of linear cones2 
    *l = *l + num_int;
}


/* 
 * Updates the solver's lb and ub contraints for integer variables
 * to the bounds specified by the node
 */
void modify_constraints( misocp_pwork* prob, misocp_node* node){
    memcpy( node->constraints, prob->ecos_prob->b, 
        prob->num_int *2*sizeof(pfloat));
}

int main()
{
    
    return 0;
}
