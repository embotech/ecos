#include "ecos.h"
#include "misocp.h"
#include "math.h"
#include "stdlib.h"
#include "splamm.h"

/* Augments the G and b arrays to take lb and ub constraints
 * for all of the variables marked integer
 * USES MALLOC
 */ 
void socp_to_misocp(
    idxint num_int, idxint n, idxint m,
    pfloat* Gpr_in, idxint* Gjc_in, idxint* Gir_in,
    pfloat* Gpr_out, idxint* Gjc_out, idxint* Gir_out,
    pfloat* h_in, pfloat* h_out)
{
    idxint i, j;
    
    spmat* G_old = createSparseMatrix(m, n, Gjc_in[n], 
        Gjc_in, Gir_in, Gpr_in);
    printSparseMatrix(G_old);
    printf("======\n\n");

    // Fill in column pointers making room for the lb's and ub's
    for (i=1; i<=n ; ++i){
        idxint a = i < num_int ? i : num_int;
        Gjc_out[i] = Gjc_in[i] + 2 * a; 
    }

    spmat* G_new = createSparseMatrix(m + 2*num_int, n, Gjc_out[n], 
        Gjc_out, Gir_out, Gpr_out);
    printSparseMatrix(G_new);    
    printf("======\n\n");

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
    
            // Set lower bound to -Inf
            h_out[ 2*i ] = INFINITY;     

            // Set upper bound to +Inf
            h_out[ 2*i + 1] = INFINITY;     

        } else {
            for(j=0; j<(Gjc_in[i+1] - Gjc_in[i]); ++j){
                Gpr_out[Gjc_out[i]+j] = Gpr_in[Gjc_in[i]+j];
                Gir_out[Gjc_out[i]+j] = Gir_in[Gjc_in[i]+j] + 2*num_int;
            }
        }
    }

    printSparseMatrix(G_new);    
    printf("======\n\n");

    // Copy the remaining entries of b
    for (i=0; i<m; ++i){
        h_out[2*num_int+i] = h_in[i];
    }
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
   
    pfloat* Gpr_new, * h_new;
    idxint* Gjc_new, * Gir_new; 

    idxint new_G_size = Gjc[n] + 2 * num_int_vars;
    Gpr_new = (pfloat *) malloc( new_G_size * sizeof(pfloat) );
    Gjc_new = (idxint *) malloc( (n+1) * sizeof(idxint) );
    Gir_new = (idxint *) malloc( new_G_size * sizeof(idxint) );
    h_new = (pfloat *) malloc( (m + 2 * num_int_vars ) * sizeof(pfloat) );

    // Copy the data and convert it to boolean
    socp_to_misocp(num_int_vars, n, m,
                    Gpr, Gjc, Gir,
                    Gpr_new, Gjc_new, Gir_new,
                    h, h_new);
    m += 2*num_int_vars;
    l += 2*num_int_vars;

    for (i=0; i<m; ++i){
        printf("%f,", h_new[l]);
    }

    // Malloc the problem's memory
    misocp_pwork* prob = (misocp_pwork*) malloc(sizeof(misocp_pwork));

    // Malloc the best optimal solution's memory
    prob->best_x = (pfloat*) malloc( n*sizeof(pfloat) );

    // Malloc the initial node's book keeping data #(2*maxIter)
    prob->nodes = (node*) calloc( 2*MI_MAXITER, sizeof(node) );

    // Malloc the id array
    prob->node_ids = (char*) malloc( 2*MI_MAXITER*num_int_vars*sizeof(char) );

    printf("\nPassed mallocs \n");

    // Setup the ecos solver
    prob->ecos_prob = ECOS_setup(
        n, m, p, l, ncones, q,
        Gpr_new, Gjc_new, Gir_new,
        Apr, Ajc, Air,
        c, h_new, b);

    printf("\nPassed ECOS \n");

    // Store the number of integer variables in the problem
    prob->num_int_vars = num_int_vars;

    prob->global_U = INFINITY;

    // Now initialize the root node
    prob->nodes[0].status = MI_NOT_SOLVED;
    prob->nodes[0].L = -INFINITY;
    prob->nodes[0].U =  INFINITY;
    for (i=0; i<num_int_vars; ++i){ prob->node_ids[i] = MI_STAR; }

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



