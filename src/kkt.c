/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
 * Automatic Control Laboratory, ETH Zurich.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/* The KKT module.
 * Handles all computation related to KKT matrix:
 * - updating the matrix
 * - its factorization
 * - solving for search directions
 * - etc.
 */

#include "kkt.h"
#include "ldl.h"
#include "splamm.h"
#include "ecos.h"
#include "cone.h"

#include <math.h>

/* Factorization of KKT matrix. Just a wrapper for some LDL code */
idxint kkt_factor(kkt* KKT, pfloat eps, pfloat delta)
{
	idxint nd;
    
    /* returns n if successful, k if D (k,k) is zero */
	nd = LDL_numeric2(
				KKT->PKPt->n,	/* K and L are n-by-n, where n >= 0 */
				KKT->PKPt->jc,	/* input of size n+1, not modified */
				KKT->PKPt->ir,	/* input of size nz=Kjc[n], not modified */
				KKT->PKPt->pr,	/* input of size nz=Kjc[n], not modified */
				KKT->L->jc,		/* input of size n+1, not modified */
				KKT->Parent,	/* input of size n, not modified */
				KKT->Sign,      /* input, permuted sign vector for regularization */
                eps,            /* input, inverse permutation vector */
				delta,          /* size of dynamic regularization */
				KKT->Lnz,		/* output of size n, not defn. on input */
				KKT->L->ir,		/* output of size lnz=Lp[n], not defined on input */
				KKT->L->pr,		/* output of size lnz=Lp[n], not defined on input */
				KKT->D,			/* output of size n, not defined on input */
				KKT->work1,		/* workspace of size n, not defn. on input or output */
				KKT->Pattern,   /* workspace of size n, not defn. on input or output */
				KKT->Flag	    /* workspace of size n, not defn. on input or output */				
    );

	return nd == KKT->PKPt->n ? KKT_OK : KKT_PROBLEM;
}


/**
 * Solves the permuted KKT system and returns the unpermuted search directions.
 *
 * On entry, the factorization of the permuted KKT matrix, PKPt,
 * is assumed to be up to date (call kkt_factor beforehand to achieve this).
 * The right hand side, Pb, is assumed to be already permuted.
 *
 * On exit, the resulting search directions are written into dx, dy and dz,
 * where these variables are permuted back to the original ordering.
 *
 * KKT->nitref iterative refinement steps are applied to solve the linear system.
 *
 * Returns the number of iterative refinement steps really taken.
 */
idxint kkt_solve(kkt* KKT, spmat* A, spmat* G, pfloat* Pb, pfloat* dx, pfloat* dy, pfloat* dz, idxint n, idxint p, idxint m, cone* C, idxint isinit, idxint nitref)
{
	idxint i, k, l, j, kk, kItRef;
	idxint*  Pinv = KKT->Pinv;
	pfloat*    Px = KKT->work1;
	pfloat*   dPx = KKT->work2;
	pfloat*     e = KKT->work3;
    pfloat*    Pe = KKT->work4;
    pfloat* truez = KKT->work5;
    pfloat*   Gdx = KKT->work6;
    pfloat* ex = e;
    pfloat* ey = e + n;
    pfloat* ez = e + n+p;
    pfloat bnorm = 1.0 + norminf(Pb, n+p+m+2*C->nsoc);
    pfloat nex, ney, nez;
    pfloat error_threshold = bnorm*LINSYSACC;
    idxint nK = KKT->PKPt->n;

	/* forward - diagonal - backward solves: Px holds solution */		
	LDL_lsolve2(nK, Pb, KKT->L->jc, KKT->L->ir, KKT->L->pr, Px );		
	LDL_dsolve(nK, Px, KKT->D);
	LDL_ltsolve(nK, Px, KKT->L->jc, KKT->L->ir, KKT->L->pr);
    
#if PRINTLEVEL > 2
    PRINTTEXT("\nIR: it  ||ex||   ||ey||   ||ez|| (threshold: %4.2e\n", error_threshold);
    PRINTTEXT("    -------------------------------------------------\n");
#endif
    
	/* iterative refinement */
	for( kItRef=1; kItRef <= nitref; kItRef++ ){
        
        /* unpermute x & copy into arrays */
        unstretch(n, p, m, C, Pinv, Px, dx, dy, dz);
        
		/* compute error term */
        k=0; j=0;
        /* --> 1. ex = bx - A'*dy - G'*dz */
        for( i=0; i<n; i++ ){ ex[i] = Pb[Pinv[k++]]; }
        sparseMtVm(A, dy, ex, 0, 0);
        sparseMtVm(G, dz, ex, 0, 0);
        nex = norminf(ex,n);
        	
        /* --> 2. ey = by - A*dx */
        for( i=0; i<p; i++ ){ ey[i] = Pb[Pinv[k++]]; }
        sparseMV(A, dx, ey, -1, 0);
        ney = norminf(ey,p);
        
        /* --> 3. ez = bz - G*dx + V*dz_true */
        kk = 0;
        sparseMV(G, dx, Gdx, 1, 1);
        j=0;
        for( i=0; i<C->lpc->p; i++ ){
            ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
        }
        for( l=0; l<C->nsoc; l++ ){
            for( i=0; i<C->soc[l].p; i++ ){
                ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
            }
            ez[kk] = 0;
            ez[kk+1] = 0;
            k += 2;
            kk += 2;
        }
        for( i=0; i<m+2*C->nsoc; i++) { truez[i] = Px[Pinv[n+p+i]]; }
        if( isinit == 0 ){
            scale2add(truez, ez, C);
        } else {
            vadd(m+2*C->nsoc, truez, ez);
        }
        nez = norminf(ez,m+2*C->nsoc);
        
        
#if PRINTLEVEL > 2
        PRINTTEXT("    %2d  %3.1e  %3.1e  %3.1e\n", (int)kItRef, nex/bnorm, ney/bnorm, nez/bnorm);
#endif
        
        /* continue with refinement only if errors are small enough */
        if( nex < error_threshold && ney < error_threshold && nez < error_threshold){
            kItRef--;
            break;
        }
        
        /* permute */
        for( i=0; i<nK; i++) { Pe[Pinv[i]] = e[i]; }
        
        /* forward - diagonal - backward solves: dPx holds solution */
        LDL_lsolve2(nK, Pe, KKT->L->jc, KKT->L->ir, KKT->L->pr, dPx);
        LDL_dsolve(nK, dPx, KKT->D);
        LDL_ltsolve(nK, dPx, KKT->L->jc, KKT->L->ir, KKT->L->pr);
        
        /* add refinement to Px */
        for( i=0; i<nK; i++ ){ Px[i] += dPx[i]; }
	}

#if PRINTLEVEL > 2
    PRINTTEXT("\n");
#endif
    
	/* copy solution out into the different arrays, permutation included */
	unstretch(n, p, m, C, Pinv, Px, dx, dy, dz);
    
    return kItRef == nitref+1 ? nitref : kItRef;
}


/**
 * Updates the permuted KKT matrix by copying in the new scalings.
 */
void kkt_update(spmat* PKP, idxint* P, cone *C)
{
	idxint i, j, k, conesize, conesize_m1;
    pfloat eta_square, d1, u0, u1, v1, *q;
	
	/* LP cone */
    for( i=0; i < C->lpc->p; i++ ){ PKP->pr[P[C->lpc->kkt_idx[i]]] = -C->lpc->v[i]; }

	/* Second-order cone */
	for( i=0; i<C->nsoc; i++ ){
        
        getSOCDetails(&C->soc[i], &conesize, &eta_square, &d1, &u0, &u1, &v1, &q);
        conesize_m1 = conesize - 1;
        
        /* D */
        PKP->pr[P[C->soc[i].Didx[0]]] = -eta_square * d1;
        for (k=1; k < conesize; k++) {
            PKP->pr[P[C->soc[i].Didx[k]]] = -eta_square;
        }
        
        
        /* v */
        j=1;
        for (k=0; k < conesize_m1; k++) {
            PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square * v1 * q[k];
        }
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square;
        
        /* u */
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square * u0;
        for (k=0; k < conesize_m1; k++) {
            PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square * u1 * q[k];
        }
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = +eta_square;
	}
}
