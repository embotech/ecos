#include "ecos.h"
#include "misocp.h"
#include "math.h"
#include "stdlib.h"

// Same as string.h's memcpy
void memory_copy(void* src, void* dest, idxint bytes){
    int i;
    for (i=0; i<bytes; ++i){
        ((char*) dest)[i] = ((char*) src)[i];
    }
}

/* Augments the G and b arrays to take lb and ub constraints
 * for all of the variables marked integer
 * USES MALLOC
 */ 
void socp_to_misocp(
    idxint num_int, idxint n, idxint m, idxint* l,
    pfloat* Gpr_in, idxint* Gjc_in, idxint* Gir_in,
    pfloat* Gpr_out, idxint* Gjc_out, idxint* Gir_out,
    pfloat* b_in, pfloat* b_out)
{
    idxint i;
    /** We assume that the first num_int column vars are integers **/
    
    // First we need to resize the G mat to accomodate lb and ub constraints
    idxint new_G_size = Gjc_in[n] + 2 * num_int;
    Gpr_out = (pfloat *) malloc( new_G_size * sizeof(pfloat) );
    Gjc_out = (idxint *) malloc( n * sizeof(idxint) );
    Gir_out = (idxint *) malloc( new_G_size * sizeof(idxint) );

    b_out = (pfloat *) malloc( (n + 2 * num_int ) * sizeof(pfloat) );

    // Fill in column pointers making room for the lb's and ub's
    for (i=0; i<n ; ++i){
        idxint a = i < num_int ? i : num_int;
        Gjc_out[i] = Gjc_in[i] + 2 * a; 
    }

    // Fill in the new matrices 
    for (i=0; i<num_int; ++i){
        if (i<num_int){
            // Set coefficients for G matrix
            Gpr_out[ Gjc_out[i] ] = -1;
            Gpr_out[ Gjc_out[i] + 1 ] = 1;

            Gir_out[ Gjc_out[i] ] = 2*i;
            Gir_out[ Gjc_out[i] + 1 ] = 2*i + 1;

            // Copy data
            memory_copy(&Gpr_in[ Gjc_in[i] ] , &Gpr_out[ Gjc_out[i] + 2],
                ( Gjc_in[i+1] - Gjc_in[i] )*sizeof(pfloat) );        

            // Copy indexes
            memory_copy(&Gir_in[ Gjc_in[i] ] , &Gir_out[ Gjc_out[i] + 2],
                ( Gjc_in[i+1] - Gjc_in[i] )*sizeof(idxint) );        

            // Set upper bound to +Inf
            b_out[ 2*i ] = INFINITY;     

            // Set lower bound to -Inf
            b_out[ 2*i + 1] = -INFINITY;     
        } else {
            // Copy data
            memory_copy(&Gpr_in[ Gjc_in[i] ] , &Gpr_out[ Gjc_out[i] ],
                ( Gjc_in[i+1] - Gjc_in[i] )*sizeof(pfloat) );        

            // Copy indexes
            memory_copy(&Gir_in[ Gjc_in[i] ] , &Gir_out[ Gjc_out[i] ],
                ( Gjc_in[i+1] - Gjc_in[i] )*sizeof(idxint) );        

        }
    }

    // Copy the remaining entries of b
    memory_copy( b_in, &b_out[2*num_int], m*sizeof(pfloat) );

    // Augment the number of linear cones2 
    *l = *l + num_int;
}


// We assume that the first num_int column vars are booleans
// Boolean vars only
misocp_pwork* misocp_setup(
    idxint n, idxint m, idxint p, 
    idxint l, idxint ncones, idxint* q,
    pfloat* Gpr, idxint* Gjc, idxint* Gir,
    pfloat* Apr, idxint* Ajc, idxint* Air,
    pfloat* c, pfloat* h, pfloat* b, idxint num_int_vars)
{
    idxint i;
   
    pfloat* Gpr_new, * b_new;
    idxint* Gjc_new, * Gir_new; 

    // Copy the data and convert it to boolean
    socp_to_misocp(num_int_vars, n, m, &l,
                    Gpr, Gjc, Gir,
                    Gpr_new, Gjc_new, Gir_new,
                    b, b_new);

    //idxint num_int, idxint n, idxint m, idxint* l,
    //pfloat* Gpr_in, idxint* Gjc_in, idxint* Gir_in,
    //pfloat* Gpr_out, idxint* Gjc_out, idxint* Gir_out,
    //pfloat* b_in, pfloat* b_out)

    // Malloc the problem's memory
    misocp_pwork* prob = (misocp_pwork*) malloc(sizeof(misocp_pwork));

    // Setup the ecos solver
    prob->ecos_prob = ECOS_setup(
        n, m, p, l, ncones, q,
        Gpr, Gjc, Gir,
        Apr, Ajc, Air,
        c, h, b);

    // Malloc the best optimal solution's memory
    prob->best_x = malloc( n*sizeof(pfloat) );

    // Malloc the initial node's book keeping data #(2*maxIter)
    prob->nodes = (node*) calloc( 2*MI_MAXITER, sizeof(node) );

    // Malloc the id array
    prob->node_ids = (char*) malloc( 2*MI_MAXITER*num_int_vars*sizeof(char) );

    // Store the number of integer variables in the problem
    prob->num_int_vars = num_int_vars;

    prob->global_U = INFINITY;

    // Now initialize the root node
    prob->nodes[0].status = MI_NOT_SOLVED;
    prob->nodes[0].L = -INFINITY;
    prob->nodes[0].U =  INFINITY;
    for (i=0; i<num_int_vars; ++i){
        prob->node_ids[i] = MI_STAR;
    }

    return prob;
}



// Performs the same function as ecos_cleanup
void misocp_cleanup(misocp_pwork* prob){
    // Free solver memory
    ECOS_cleanup(prob->ecos_prob, 0);

    free(prob->nodes);
    free(prob->node_ids);
    free(prob->best_x);
    free(prob);
}



