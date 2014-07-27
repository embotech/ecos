#include "ecos.h"
#include "ecos_bb.h"
#include "math.h"
#include "stdlib.h"
#include "splamm.h"

/* Augments the G and b arrays to take lb and ub constraints
 * for all of the variables marked integer
 * USES MALLOC
 */ 
void socp_to_ecos_bb(
    idxint num_int, idxint n, idxint m,
    pfloat* Gpr_in, idxint* Gjc_in, idxint* Gir_in,
    pfloat* Gpr_out, idxint* Gjc_out, idxint* Gir_out,
    pfloat* h_in, pfloat* h_out)
{
    idxint i, j;

    // Fill in column pointers making room for the lb's and ub's
    Gjc_out[0] = 0;
    for (i=1; i<=n ; ++i){
        idxint a = i < num_int ? i : num_int;
        Gjc_out[i] = Gjc_in[i] + 2 * a; 
    }

    // Fill in the new matrices 
    for (i=0; i<n; ++i){
        if (i<num_int){
            // Set coefficients for G matrix
            Gpr_out[ Gjc_out[i] ] = -1;
            Gpr_out[ Gjc_out[i] + 1 ] = 1;

            Gir_out[ Gjc_out[i] ] = 2*i;
            Gir_out[ Gjc_out[i] + 1 ] = 2*i + 1;

            for(j=0; j<(Gjc_in[i+1] - Gjc_in[i]); ++j){
                Gpr_out[Gjc_out[i]+2+j] = Gpr_in[Gjc_in[i]+j];
                Gir_out[Gjc_out[i]+2+j] = Gir_in[Gjc_in[i]+j] + 2*num_int;
            }
            
            // Set lower bound to 0
            h_out[ 2*i ] = 0;     

            // Set upper bound to 1
            h_out[ 2*i + 1] = 1;     

        } else {
            for(j=0; j<(Gjc_in[i+1] - Gjc_in[i]); ++j){
                Gpr_out[Gjc_out[i]+j] = Gpr_in[Gjc_in[i]+j];
                Gir_out[Gjc_out[i]+j] = Gir_in[Gjc_in[i]+j] + 2*num_int;
            }
        }
    }

    // Copy the remaining entries of b
    for (i=0; i<m; ++i){
        h_out[2*num_int+i] = h_in[i];
    }
}


/* We assume that the first num_int column vars are booleans, otherwise all 
 * arguements are exactly the same as ECOS
 * Boolean vars only
*/
ecos_bb_pwork* ecos_bb_setup(
    idxint n, idxint m, idxint p, 
    idxint l, idxint ncones, idxint* q,
    pfloat* Gpr, idxint* Gjc, idxint* Gir,
    pfloat* Apr, idxint* Ajc, idxint* Air,
    pfloat* c, pfloat* h, pfloat* b, idxint num_bool_vars)
{
    
    pfloat* Gpr_new, * h_new;
    idxint* Gjc_new, * Gir_new; 

    idxint new_G_size = Gjc[n] + 2 * num_bool_vars;
    Gpr_new = (pfloat *) malloc( new_G_size * sizeof(pfloat) );
    Gjc_new = (idxint *) malloc( (n+1) * sizeof(idxint) );
    Gir_new = (idxint *) malloc( new_G_size * sizeof(idxint) );
    h_new = (pfloat *) malloc( (m + 2 * num_bool_vars ) * sizeof(pfloat) );

    // Copy the data and convert it to boolean
    socp_to_ecos_bb(num_bool_vars, n, m,
                    Gpr, Gjc, Gir,
                    Gpr_new, Gjc_new, Gir_new,
                    h, h_new);
    m += 2*num_bool_vars;
    l += 2*num_bool_vars;

    // Malloc the problem's memory
    ecos_bb_pwork* prob = (ecos_bb_pwork*) malloc(sizeof(ecos_bb_pwork));

    // Malloc the best optimal solution's memory
    prob->best_x = (pfloat*) malloc( n*sizeof(pfloat) );

    // Malloc the initial node's book keeping data #(2*maxIter)
    prob->nodes = (node*) calloc( MI_MAXITER, sizeof(node) );

    // Malloc the id array
    prob->node_ids = (char*) malloc( MI_MAXITER*num_bool_vars*sizeof(char) );

    // Malloc the tmp node id
    prob->tmp_node_id = (char*) malloc( num_bool_vars*sizeof(char) );

    // Setup the ecos solver
    prob->ecos_prob = ECOS_setup(
        n, m, p, l, ncones, q,
        Gpr_new, Gjc_new, Gir_new,
        Apr, Ajc, Air,
        c, h_new, b);

    // Store the number of integer variables in the problem
    prob->num_bool_vars = num_bool_vars;

    prob->global_U = INFINITY;

    // offset the h pointer for the user
    prob->h = &prob->ecos_prob->h[2 * num_bool_vars];

    // Map the other variables
    prob->A = prob->ecos_prob->A;
    prob->G = prob->ecos_prob->G;
    prob->c = prob->ecos_prob->c;
    prob->b = prob->ecos_prob->b;
    
    /* switch off ecos prints */
    prob->ecos_prob->stgs->verbose = 0;
    
    return prob;
}

// Performs the same function as ecos_cleanup
void ecos_bb_cleanup(ecos_bb_pwork* prob){
    // Free solver memory
    ECOS_cleanup(prob->ecos_prob, 0);
    free(prob->tmp_node_id);
    free(prob->nodes);
    free(prob->node_ids);
    free(prob->best_x);
    free(prob);
}



