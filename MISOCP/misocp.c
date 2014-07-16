#include "ecos.h"
#include "stdlib.h"
#include "math.h"
#include "misocp.h"

/* 
 * Modify the input parameters before creating the template ECOS problem
 */
misocp_pwork* miscop_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones, 
    idxint* q,
    pfloat* Gpr, idxint* Gjc, idxint* Gir,
    pfloat* Apr, idxint* Ajc, idxint* Air,
    pfloat* c, pfloat* h, pfloat* b, idxint num_int)
{
    idxint i, j;
    misocp_pwork* misocp_prob = (misocp_pwork*) malloc(sizeof(misocp_pwork));

    /** We assume that the first num_int column vars are integers **/
    
    // First we need to resize the G mat to accomodate lb and ub constraints
    idxint* new_G_size = Gjc[n] + 2 * num_int;
    pfloat* Gpr_new = (pfloat *) malloc( new_G_size * sizeof(pfloat) );
    idxint* Gjr_new = (idxint *) malloc( n * sizeof(idxint) );
    idxint* Gir_new = (idxint *) malloc( new_G_size * sizeof(idxint) );

    pfloat* b_new = (pfloat *) malloc( (n + 2 * num_int ) * sizeof(pfloat) );

    // Fill in column pointers making room for the lb's and ub's
    for (i=0; i<n ; ++i){
        Gjc_new[i] = Gjc[i] + 2 * min(i, num_int); 
    }

    // Fill in the new matrices 
    for (i=0; i<num_int; ++i){
        if (i<num_int){
            // Set coefficients for G matrix
            Gpr_new[ Gjc_new[i] ] = 1;
            Gpr_new[ Gjc_new[i] + 1 ] = -1;

            Gir_new[ Gjc_new[i] ] = 2*i;
            Gir_new[ Gjc_new[i] + 1 ] = 2*i + 1;

            // Copy data
            memcpy(&Gpr[Gjr[i]] , &Gpr_new[Gjr_new[i] + 2]
                (Gjr[i+1] - Gjr[i]) * sizeof(pfloat) );        

            // Copy indexes
            memcpy(&Gir[Gjr[i]] , &Gir_new[Gjr_new[i] + 2]
                (Gjr[i+1] - Gjr[i]) * sizeof(idxint) );        

            // Set upper bound to +Inf
            b_new[ 2*i ] = 1e999;     

            // Set lower bound to -Inf
            b_new[ 2*i + 1] = -1e999;     
        } else {
            // Copy data
            memcpy(&Gpr[Gjr[i]] , &Gpr_new[ Gjr_new[i] ]
                (Gjr[i+1] - Gjr[i]) * sizeof(pfloat) );        

            // Copy indexes
            memcpy(&Gir[Gjr[i]] , &Gir_new[ Gjr_new[i] ]
                (Gjr[i+1] - Gjr[i]) * sizeof(idxint) );        

        }
    }

    // Increase the count for linear cones
    l = l + num_int;

    // Setup the ecos solver
    misocp_prob->ecos_prob = ECOS_setup(n, m, p, l, ncones, q,
        Gpr_new, Gjc_new, Gir_new,
        Apr, Ajc, Air,
        c, h, b_new);

    return misocp_prob;
}

void modify_socp_prob( misocp_pwork* prob, misocp_node* node){
    memcpy( node->constraints, prob->ecos_prob->b, 
        prob->num_int *2*sizeof(pfloat));
}

int main()
{
    pwork *problem;
    idxint exitflag;

    /*
    pwork* ECOS_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones, idxint* q,
    pfloat* Gpr, idxint* Gjc, idxint* Gir,
    pfloat* Apr, idxint* Ajc, idxint* Air,
    pfloat* c, pfloat* h, pfloat* b);
    */

    /* set up data */
    mywork = ECOS_setup(1, 2, 0,
        2, 0, NULL,
        feas_Gx, feas_Gp, feas_Gi,
        NULL, NULL, NULL,
        feas_c, feas_h, NULL);

    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork); 
    }
    else exitflag = ECOS_FATAL;
        
    /* clean up memory */
    ECOS_cleanup(mywork, 0);
    
    mu_assert("feas-test: ECOS failed to produce outputflag OPTIMAL", exitflag == ECOS_OPTIMAL );
    return 0;
}
