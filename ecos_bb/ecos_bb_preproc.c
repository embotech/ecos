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

int contains(idxint idx, idxint num_int, idxint* bool_vars_idx){
    idxint i;
    idxint ans = 0;
    for (i=0; i<num_int; ++i){
        ans += (bool_vars_idx[i] == idx);
    }
    return ans;
}

/* Augments the G and b arrays to take lb and ub constraints
 * for all of the variables marked integer
 * USES MALLOC
 */ 
void socp_to_ecos_bb(
    idxint num_int, idxint* bool_vars_idx, 
    idxint n, idxint m,
    pfloat* Gpr_in, idxint* Gjc_in, idxint* Gir_in,
    pfloat* Gpr_out, idxint* Gjc_out, idxint* Gir_out,
    pfloat* h_in, pfloat* h_out)
{
    idxint i, j, k;

    /* First map in the column pointers, remember Gir_out[0]=0*/
    for (i=0; i<=n; ++i){
        Gjc_out[i] = Gjc_in[i];
    }

    /* Now insert the new zeros and expand the column indices as needed*/
    for (j=0; j<num_int; ++j){
        k = bool_vars_idx[j];
        Gpr_out[ Gjc_out[k] ] = -1;
        Gpr_out[ Gjc_out[k] + 1 ] = 1;

        Gir_out[ Gjc_out[k] ] = 2*j;
        Gir_out[ Gjc_out[k] + 1 ] = 2*j + 1;

        /* Set lower bound to 0*/
        h_out[ 2*j ] = 0;     

        /* Set upper bound to 1*/
        h_out[ 2*j + 1] = 1;     

        for (i=(k+1); i<=n; ++i){
            Gjc_out[i] += 2;
        }
    }

    /* Now fill in the remainder of the array*/
    for (i=0; i<n; ++i){
        if ( contains(i, num_int, bool_vars_idx) ){
            for(j=0; j<(Gjc_in[i+1] - Gjc_in[i]); ++j){
                Gpr_out[Gjc_out[i]+2+j] = Gpr_in[Gjc_in[i]+j];
                Gir_out[Gjc_out[i]+2+j] = Gir_in[Gjc_in[i]+j] + 2*num_int;
            }

        } else {
            for(j=0; j<(Gjc_in[i+1] - Gjc_in[i]); ++j){
                Gpr_out[Gjc_out[i]+j] = Gpr_in[Gjc_in[i]+j];
                Gir_out[Gjc_out[i]+j] = Gir_in[Gjc_in[i]+j] + 2*num_int;
            }
        }
    }

    /* Copy the remaining entries of b*/
    for (i=0; i<m; ++i){
        h_out[2*num_int+i] = h_in[i];
    }

#if PRINTLEVEL >= 3
    PRINTTEXT("Gjc_out: ");
    for (i=0; i<=n; ++i){
        PRINTTEXT("%u ", Gjc_out[i] );
    }
    PRINTTEXT("\n");

    PRINTTEXT("Gpr_out: ");
    for (i=0; i<Gjc_out[n]; ++i){
        PRINTTEXT("%f.2 ", Gpr_out[i] );
    }
    PRINTTEXT("\n");

    PRINTTEXT("Gir_out: ");
    for (i=0; i<Gjc_out[n]; ++i){
        PRINTTEXT("%u ", Gir_out[i] );
    }
    PRINTTEXT("\n");
#endif

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
    pfloat* c, pfloat* h, pfloat* b, 
    idxint num_bool_vars, idxint* bool_vars_idx)
{
    
    pfloat* Gpr_new, * h_new;
    idxint* Gjc_new, * Gir_new; 

    idxint new_G_size = Gjc[n] + 2 * num_bool_vars;
    Gpr_new = (pfloat *) MALLOC( new_G_size * sizeof(pfloat) );
    Gjc_new = (idxint *) MALLOC( (n+1) * sizeof(idxint) );
    Gir_new = (idxint *) MALLOC( new_G_size * sizeof(idxint) );
    h_new = (pfloat *) MALLOC( (m + 2 * num_bool_vars ) * sizeof(pfloat) );

    /* Copy the data and convert it to boolean*/
    socp_to_ecos_bb(num_bool_vars, bool_vars_idx,
                    n, m,
                    Gpr, Gjc, Gir,
                    Gpr_new, Gjc_new, Gir_new,
                    h, h_new);
    m += 2*num_bool_vars;
    l += 2*num_bool_vars;

    /* MALLOC the problem's memory*/
    ecos_bb_pwork* prob = (ecos_bb_pwork*) MALLOC(sizeof(ecos_bb_pwork));

    /* Default the maxiter to global declared in the header file*/
    prob->maxiter = MI_MAXITER;

    /* MALLOC the best optimal solution's memory*/
    prob->best_x = (pfloat*) MALLOC( n*sizeof(pfloat) );

    /* MALLOC the initial node's book keeping data #(2*maxIter)*/
    prob->nodes = (node*) CALLOC( MI_MAXITER, sizeof(node) );

    /* MALLOC the id array*/
    prob->node_ids = (char*) MALLOC( MI_MAXITER*num_bool_vars*sizeof(char) );

    /* MALLOC the tmp node id*/
    prob->tmp_node_id = (char*) MALLOC( num_bool_vars*sizeof(char) );

    /* Store the pointer to the boolean idx*/
    prob->bool_vars_idx = bool_vars_idx;

    /* Setup the ecos solver*/
    prob->ecos_prob = ECOS_setup(
        n, m, p, l, ncones, q,
        Gpr_new, Gjc_new, Gir_new,
        Apr, Ajc, Air,
        c, h_new, b);

    /* Store the number of integer variables in the problem*/
    prob->num_bool_vars = num_bool_vars;

    prob->global_U = INFINITY;

    /* offset the h pointer for the user*/
    prob->h = &prob->ecos_prob->h[2 * num_bool_vars];

    /* Map the other variables*/
    prob->A = prob->ecos_prob->A;
    prob->G = prob->ecos_prob->G;
    prob->c = prob->ecos_prob->c;
    prob->b = prob->ecos_prob->b;
    
    /* switch off ecos prints */
    prob->ecos_prob->stgs->verbose = 0;
    
    return prob;
}

/* Performs the same function as ecos_cleanup*/
void ecos_bb_cleanup(ecos_bb_pwork* prob, idxint num_vars_keep){
    /* FREE solver memory*/
    ECOS_cleanup(prob->ecos_prob, num_vars_keep);
    FREE(prob->tmp_node_id);
    FREE(prob->nodes);
    FREE(prob->node_ids);
    FREE(prob->best_x);
    FREE(prob);
}


