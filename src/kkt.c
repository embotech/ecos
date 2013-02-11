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
idxint kkt_factor(kkt* KKT, pfloat delta)
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
                KKT->Pinv,      /* input, inverse permutation vector */
				delta,			/* size of dynamic regularization */
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
 */
void kkt_solve(kkt* KKT, spmat* A, spmat* G, pfloat* Pb, pfloat* dx, pfloat* dy, pfloat* dz, idxint n, idxint p, idxint m, cone* C, idxint isinit, idxint nitref)
{
	idxint i, k, l, j, kk, nK, kItRef;
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
    pfloat Pxnorm, temp;
	nK = KKT->PKPt->n;
    pfloat bnorm = norm2(Pb, n+p+m+C->nsoc);
    pfloat nex, ney, nez;


	/* forward - diagonal - backward solves: Px holds solution */		
	LDL_lsolve2(nK, Pb, KKT->L->jc, KKT->L->ir, KKT->L->pr, Px );		
	LDL_dsolve(nK, Px, KKT->D);
	LDL_ltsolve(nK, Px, KKT->L->jc, KKT->L->ir, KKT->L->pr);
    
	/* iterative refinement */
	for( kItRef=1; kItRef <= nitref; kItRef++ ){
        
        /* unpermute x & copy into arrays */
        k = 0; j=0;
        for( i=0; i<n; i++ ){ dx[i] = Px[Pinv[k++]]; }
        for( i=0; i<p; i++ ){ dy[i] = Px[Pinv[k++]]; }
        for( i=0; i<C->lpc->p; i++ ){ dz[j++] = Px[Pinv[k++]]; }
        for( l=0; l<C->nsoc; l++ ){
            for( i=0; i<C->soc[l].p; i++ ){ dz[j++] = Px[Pinv[k++]]; }
            k++;
        }
        
		/* compute error term */
        k=0; j=0;
        /* --> 1. ex = bx - A'*dy - G'*dz */
        for( i=0; i<n; i++ ){ ex[i] = Pb[Pinv[k++]]; }
        sparseMtVm(A, dy, ex, 0, 0);
        sparseMtVm(G, dz, ex, 0, 0);
        nex = norm2(ex,n);
        //PRINTTEXT("norm(ex) = %8.6e (k=%d)\n", norm2(ex, 48), kItRef);
        //for( i=0; i<n; i++){ PRINTTEXT("ex[%d] = %6.4e\n", i+1, ex[i]); }
        	
        /* --> 2. ey = by - A*dx */
        for( i=0; i<p; i++ ){ ey[i] = Pb[Pinv[k++]]; }
        sparseMV(A, dx, ey, -1, 0);
        ney = norm2(ey,p);
        //PRINTTEXT("norm(ey) = %8.6e (k=%d)\n", norm2(ey, 35), kItRef);
        //for( i=0; i<p; i++){ PRINTTEXT("ey[%d] = %6.4e\n", i+1, ey[i]); }
        
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
            ez[kk++] =  Pb[Pinv[k++]];
        }
        for( i=0; i<m+C->nsoc; i++) { truez[i] = Px[Pinv[n+p+i]]; }
        
        
        //PRINTTEXT("dz =   \n");
        //for( i=0; i<m+C->nsoc; i++){ PRINTTEXT("%6.4e  ", truez[i]); }
        //PRINTTEXT("\n");
        
        //PRINTTEXT("b - G*dx =   \n");
        //for( i=0; i<m+C->nsoc; i++){ PRINTTEXT("%6.4e  ", ez[i]); }
        //PRINTTEXT("\n");
        
        
        if( isinit == 0 ){
            scale2add(truez, ez, C);
            //scale2add(dz, ez, C);
        } else {
            vadd(m+C->nsoc, truez, ez);
            //vadd(m+C->nsoc, dz, ez);
        }
        //PRINTTEXT("norm(ez) = %8.6e (k=%d)\n", norm2(ez, 35), kItRef);
        //if( norm2(ez,35) < DELTA ) break;
        //if( norm2(ez,35) > 1e-5 ) kItRef--;
        //PRINTTEXT("ez =   \n");
        //for( i=0; i<m+C->nsoc; i++){ PRINTTEXT("%6.4e  ", ez[i]); }
        //PRINTTEXT("\n");
        
        //if( isinit == 0 )
        //exit(-1);
        
        nez = norm2(ez,m+C->nsoc);
        
        
        /* break here if errors are small enough */
        if( nez/bnorm < LINSYSACC && ney/bnorm < LINSYSACC && nex/bnorm < LINSYSACC){
            break;
        }
        
        PRINTTEXT("*");
        
        /* permute */
        for( i=0; i<nK; i++) { Pe[Pinv[i]] = e[i]; }
        
        /* forward - diagonal - backward solves: dPx holds solution */
        LDL_lsolve2(nK, Pe, KKT->L->jc, KKT->L->ir, KKT->L->pr, dPx);
        LDL_dsolve(nK, dPx, KKT->D);
        LDL_ltsolve(nK, dPx, KKT->L->jc, KKT->L->ir, KKT->L->pr);

        //for( i=0; i<nK; i++){ PRINTTEXT("dPx[%d] = %6.4e\n", i+1, dPx[Pinv[i]]); }
        //exit(-1);
        
        /* add refinement to Px */
        for( i=0; i<nK; i++ ){ Px[i] += dPx[i]; }
	}
    PRINTTEXT("\n");
	
	/* copy solution out into the different arrays, permutation included */
	k = 0; j=0;
	for( i=0; i<n; i++ ){ dx[i] = Px[Pinv[k++]]; }
	for( i=0; i<p; i++ ){ dy[i] = Px[Pinv[k++]]; }
	for( i=0; i<C->lpc->p; i++ ){ dz[j++] = Px[Pinv[k++]]; }
	for( l=0; l<C->nsoc; l++ ){
		for( i=0; i<C->soc[l].p; i++ ){ dz[j++] = Px[Pinv[k++]]; }
		k++;
	}	
}


/**
 * Updates the permuted KKT matrix by copying in the new scalings.
 */
void kkt_update(spmat* PKP, idxint* P, cone *C, idxint* Sign)
{
	idxint i, k, size_minus_1;
    pfloat eta_times_delta;
    pfloat temp, Dii;
    pfloat cumsign;
	
	/* LP cone */
    for( i=0; i < C->lpc->p; i++ ){ PKP->pr[P[C->lpc->kkt_idx[i]]] = -C->lpc->v[i]; }

	/* Second-order cone */
	for( i=0; i<C->nsoc; i++ ){
		size_minus_1 = C->soc[i].p - 1;
        eta_times_delta = DELTA;

		/* atilde */
        PKP->pr[P[C->soc[i].kkt_atilde]] = -C->soc[i].atilde - eta_times_delta; // regularized

		/* q' in U and -eta on diagonal */
        temp = -C->soc[i].qtildefact*C->soc[i].qtildefact/C->soc[i].atilde;
        cumsign = 0;
		for( k=0; k < size_minus_1; k++ ){
            
            /* precompute sign */
            Dii = 1.0 + temp*C->soc[i].qtilde[k]*C->soc[i].qtilde[k];
            temp = temp / Dii;
            cumsign += Dii < 0 ? 1.0 : 0.0;
            
            /* update scalings (with regularization) */
			PKP->pr[P[C->soc[i].kkt_qtU[k]]]   = -C->soc[i].qtilde[k];
            //PKP->pr[P[C->soc[i].kkt_qtU[k]+1]] = -C->soc[i].eta;
            PKP->pr[P[C->soc[i].kkt_qtU[k]+1]] = Dii < 0 ? -C->soc[i].eta + eta_times_delta : -C->soc[i].eta - eta_times_delta; // regularized            
		}

		/* v part in U */
		for( k=0; k < size_minus_1; k++ ){
			PKP->pr[P[C->soc[i].kkt_vU+k]] = -C->soc[i].vtilde[k];
		}

		/* -1 part on diagonal */
		PKP->pr[P[C->soc[i].kkt_vU + size_minus_1]] = C->soc[i].eta - (2*cumsign-1.0)*eta_times_delta;
	}
}
