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


/*
 * THIS FILE PRE-PROCESSES THE PROBLEM BEFORE ACTUAL SOLVE.
 * Ideally, this is the only file using any dynamic memory allocation.
 * The main tasks are to set up the data structure, initializie memory,
 * and to compute the orderings on the regularized KKT matrix.
 */
#include "ecos.h"
#include "splamm.h"

/* NEEDED FORM MEMORY ALLOCATION --------------------------------------- */
#include <stdlib.h>

/* NEEDED FOR SQRT in NORM2 -------------------------------------------- */
#include <math.h>

/* MATRIX ORDERING LIBRARY --------------------------------------------- */
#include "amd.h"
#include "amd_internal.h"

/* SPARSE LDL LIBRARY -------------------------------------------------- */
#include "ldl.h"


/* CHOOSE RIGHT MEMORY MANAGER ----------------------------------------- */
#ifdef MATLAB_MEX_FILE
#define MALLOC mxMalloc
#define FREE mxFree
#else 
#define MALLOC malloc
#define FREE free
#endif


/* PRIVATE METHODS ----------------------------------------------------- */
/* use the define below if you do not include amd_internal.h, where this
 * MAX macro is already defined. */
/* #define MAX(X, Y)  ((X) < (Y) ? (Y) : (X)) */


/*
 * Builds KKT matrix.
 * We store and operate only on the upper triangular part.
 * Replace by or use in codegen.
 *
 * We do not regularize on the multipliers/scalings.
 */
spmat* createKKT_U(spmat* Gt, spmat* At, pfloat delta, cone* C, idxint* Sign, idxint* DiagIdx)
{
	idxint i, j, k, l, r, row_stop, row, cone_strt, ks;
	idxint n = Gt->m;
	idxint m = Gt->n;
	idxint p = At->n;
	idxint nK, nnzK;
	pfloat* Kpr;
	idxint* Kjc;
	idxint* Kir;
	
	/* calculate non-zeros needed for K and the number of columns */	
	nK = n + p + m + C->nsoc;
	nnzK = n + p + m + At->nnz + Gt->nnz + C->nsoc;
	for( i=0; i<C->nsoc; i++ ){ 
		nnzK += 2*(C->soc[i].p - 1);
	}
	
	/* allocate memory for K */	
	Kpr = (pfloat *)MALLOC(nnzK*sizeof(pfloat));
	Kir = (idxint *)MALLOC(nnzK*sizeof(idxint));
	Kjc = (idxint *)MALLOC((nK+1)*sizeof(idxint));
	
	/* fill upper triangular part of K with values */	
	k = 0; ks=0;
	/* (1,1) block: regularization with +delta*I */	
	for( j=0; j<n; j++ ){
		Kjc[j] = j;
		Kir[k] = j;
		Kpr[k] = +delta;
		Sign[ks] = +1;
		DiagIdx[ks] = k;
		k++;
		ks++;
	}

	/* (1,2) and (2,2) block: [A'; -delta*I] */
	i = 0; 
	for( j=0; j<p; j++ ){
		row = At->jc[j];
		row_stop = At->jc[j+1];
		if( row <= row_stop ){
			Kjc[n+j] = k;
			while( row++ < row_stop ){
				Kir[k] = At->ir[i];
				Kpr[k++] = At->pr[i++];
			}
		}
		Kir[k] = n+j;
		Kpr[k] = -delta;
		Sign[ks] = -1;
		DiagIdx[ks] = k;
		k++;
		ks++;
	}

	/* (1,3) and (3,3) block: [G'; 0; -I] */
	/* LP cone */
	i = 0; 
	for( j=0; j < C->lpc->p; j++ ){
		row = Gt->jc[j];
		row_stop = Gt->jc[j+1];
		if( row <= row_stop ){
			Kjc[n+p+j] = k;
			while( row++ < row_stop ){
				Kir[k] = Gt->ir[i];
				Kpr[k++] = Gt->pr[i++];
			}
		}
		C->lpc->kkt_idx[j] = k;
		Kir[k] = n+p+j;		
		Kpr[k++] = -1.0;
		Sign[ks++] = 0;
	}
    
    
	/* Second-order cones */	
	cone_strt = C->lpc->p;
    for( l=0; l < C->nsoc; l++ ){
			 
		for( j=0; j <  C->soc[l].p; j++ ){
            
            /* copy in the G' copy part */
			row = Gt->jc[cone_strt+j];
			row_stop = Gt->jc[cone_strt+j+1];
			if( row <= row_stop ){
				Kjc[n+p+cone_strt+l+j] = k;
				while( row++ < row_stop ){
					Kir[k] = Gt->ir[i];
					Kpr[k++] = Gt->pr[i++];
				}
			}
            
			/* we now copy in the scaling matrix which
             * looks as follows:
             *
             *    * * * * * * * 0
             *    * *           *
             *    *   *         *         [ a  q'  0 ]      a: scalar
             *    *     *       *      =  [ q  I   v ]    q,v: vectors of size conedim- 1
             *    *       *     *         [ 0  v' -1 ]      I: identity of size conedim - 1
             *    *         *   *
             *    *           * *
             *    0 * * * * * * -1
             *
             * NOTE: only the upper triangular part (with the diagonal elements)
             *       is copied in here.
             */
            
			/* copy in q' (row above diagonal) */
			if( j>0 ){
				C->soc[l].kkt_qtU[j-1] = k;
				Kir[k] = n+p+cone_strt+l;
				Kpr[k++] = 1.0;
                
                C->soc[l].sidx[j-1] = ks;
                Sign[ks] = 0; // to be computed during kkt update
                ks++;
			} else {
                C->soc[l].kkt_atilde = k;
                Sign[ks++] = 0;//+1;
            }

			/*
             * copy in the diagonal part, i.e. a or the identity
             * we fill this -1 just to make sure there is something
             * for the symbolic factorization - true values are filled
             * in later
             */
			Kir[k] = n+p+cone_strt+l+j;
			Kpr[k++] = -1.0;
            //Sign[ks++] = 0; // to be computed during kkt update
		}
        
		/*
         * now do the last column,
		 * this is the v part (row below diagonal)
         */
		Kjc[n+p+cone_strt+l+j] = k;
		C->soc[l].kkt_vU = k;
		for( r=1; r< C->soc[l].p; r++ ){
			Kir[k] = n+p+cone_strt+l+r;
			Kpr[k++] = 1.0;
		}
        
		/* this is the lower -1 part */
		Kir[k] = n+p+cone_strt+l+r;
		Kpr[k++] = -1.0;
		Sign[ks] = 0; // to be computed during kkt update
        C->soc[l].sidx[C->soc[l].p-1] = ks;
        ks++;

        /* prepare index for next cone */
		cone_strt += C->soc[l].p;
	}
    
    PRINTTEXT("SETUP: nK=%d and ks=%d\n",nK,ks);

	/* return KKT matrix */
	return createSparseMatrix(nK, nK, nnzK, Kjc, Kir, Kpr);
}


/**
 * Cleanup: free memory (not used for embedded solvers, only standalone)
 *
 * Use the second argument to give the number of variables to NOT free.
 * This is useful if you want to use the result of the optimization without
 * copying over the arrays. One use case is the MEX interface, where we 
 * do not want to free x,y,s,z (depending on the number of LHS).
 */
void ECOS_cleanup(pwork* w, idxint keepvars)
{
	idxint i;
	
	/* Free KKT related memory */
	FREE(w->KKT->D);
	FREE(w->KKT->dx1);
	FREE(w->KKT->dx2);
	FREE(w->KKT->dy1);
	FREE(w->KKT->dy2);
	FREE(w->KKT->dz1);
	FREE(w->KKT->dz2);
	FREE(w->KKT->Flag);
	freeSparseMatrix(w->KKT->L);
	FREE(w->KKT->Lnz);
	FREE(w->KKT->Parent);
	FREE(w->KKT->Pattern);
	FREE(w->KKT->Sign);
	FREE(w->KKT->Pinv);
	FREE(w->KKT->PK);
	freeSparseMatrix(w->KKT->PKPt);
	FREE(w->KKT->RHS1);
	FREE(w->KKT->RHS2);
	FREE(w->KKT->work1);
	FREE(w->KKT->work2);
    FREE(w->KKT->work3);
    FREE(w->KKT->work4);
    FREE(w->KKT->work5);
    FREE(w->KKT->work6);
	FREE(w->KKT);

	/* Free memory for cones */
	if( w->C->lpc->p > 0 ){
		FREE(w->C->lpc->kkt_idx);
		FREE(w->C->lpc->v);
		FREE(w->C->lpc->w);
		FREE(w->C->lpc);
	}
	for( i=0; i < w->C->nsoc; i++ ){
		FREE(w->C->soc[i].kkt_qtU);
        FREE(w->C->soc[i].sidx);
		FREE(w->C->soc[i].q);
		FREE(w->C->soc[i].qtilde);
		FREE(w->C->soc[i].skbar);
		FREE(w->C->soc[i].zkbar);
		FREE(w->C->soc[i].vtilde);		
	}
	if( w->C->nsoc > 0 ){
		FREE(w->C->soc);
	}
	FREE(w->C);

	/* free stuff from pwork */
    FREE(w->W_times_dzaff);
	FREE(w->dsaff_by_W);
    FREE(w->dzaff);
    FREE(w->dsaff);
    FREE(w->zaff);
    FREE(w->saff);
	FREE(w->info);
	FREE(w->lambda);
	FREE(w->rx);
	FREE(w->ry);
	FREE(w->rz);
	FREE(w->stgs);
	FREE(w->G);	
	if( w->p > 0 ) FREE(w->A);	
	if( keepvars < 4 ) FREE(w->z);
	if( keepvars < 3 ) FREE(w->s);
	if( keepvars < 2 ) FREE(w->y);
	if( keepvars < 1 ) FREE(w->x);
	FREE(w);
}


/* 
 * Sets the (3,3) block in KKT matrix to -I by
 * overwriting the non-diagonal parts of the conic
 * multipliers with 0.
 * Assumes that the diagonal is already -I.
 */
void prepareKKT4init(spmat* K, cone* C)
{	
	idxint i, k, size_minus_1;
	for( i=0; i<C->nsoc; i++ ){
		size_minus_1 = C->soc[i].p - 1;

		/* q part in L */
		for( k=0; k < size_minus_1; k++ ){
			K->pr[C->soc[i].kkt_qL+k] = 0;
		}

		/* q' part in U and L */
		for( k=0; k < size_minus_1; k++ ){
			K->pr[C->soc[i].kkt_qtU[k]] = 0;
			K->pr[C->soc[i].kkt_qtU[k]+2] = 0;
		}

		/* q part in U */
		for( k=0; k < size_minus_1; k++ ){
			K->pr[C->soc[i].kkt_vU+k] = 0;
		}
	}
}


/* 
 * Sets the (3,3) block in KKT matrix to -I by
 * overwriting the non-diagonal parts of the conic
 * multipliers with 0.
 * Assumes that the diagonal is already -I.
 * Works only on the upper part of K.
 */
void prepareKKT4initU(spmat* K, cone* C)
{	
	idxint i, k, size_minus_1;
	for( i=0; i<C->nsoc; i++ ){
		size_minus_1 = C->soc[i].p - 1;

		/* q' part in U */
		for( k=0; k < size_minus_1; k++ ){
			K->pr[C->soc[i].kkt_qtU[k]] = 0;			
		}

		/* q part in U */
		for( k=0; k < size_minus_1; k++ ){
			K->pr[C->soc[i].kkt_vU+k] = 0;
		}
	}
}


/*
 * Sets up all data structures needed.
 * Replace by codegen
 */
pwork* ECOS_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones, idxint* q,
                   pfloat* Gpr, idxint* Gjc, idxint* Gir,
                   pfloat* Apr, idxint* Ajc, idxint* Air,
                   pfloat* c, pfloat* h, pfloat* b)
{
    idxint i, j, k, cidx, conesize, lnz, amd_result, nK, *Ljc, *Lir, *P, *Pinv, *Sign, *DiagIdx;
    pwork* mywork;
	double Control [AMD_CONTROL], Info [AMD_INFO];		
	pfloat rx, ry, rz, *Lpr;
	spmat *At, *Gt, *KU;

#if PROFILING > 0
	timer tsetup;	
#endif

#if PROFILING > 1
	timer tcreatekkt;
	timer tmattranspose;
	timer tordering;
#endif

#if PROFILING > 0
	tic(&tsetup);
#endif
   
#if PRINTLEVEL > 2
	PRINTTEXT("\n");		
	PRINTTEXT("  ****************************************************************************\n");
	PRINTTEXT("  * ECOS: Embedded Conic Solver - Sparse Interior Point method for SOCPs     *\n");
	PRINTTEXT("  *                                                                          *\n");
	PRINTTEXT("  * NOTE: The solver is heavily based on L. Vandenberghe's 'The CVXOPT       *\n");
	PRINTTEXT("  *       linear and quadratic cone program solvers', March 20, 2010.        *\n");
	PRINTTEXT("  *       [http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf]        *\n");
	PRINTTEXT("  *                                                                          *\n");
	PRINTTEXT("  *       This code uses T.A. Davis' sparse LDL package and AMD code.        *\n");
	PRINTTEXT("  *       [http://www.cise.ufl.edu/research/sparse]                          *\n");
	PRINTTEXT("  *                                                                          *\n");
	PRINTTEXT("  *       Written during a summer visit at Stanford University with S. Boyd. *\n");
	PRINTTEXT("  *                                                                          *\n");
	PRINTTEXT("  * (C) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2012.  *\n");
	PRINTTEXT("  *     Email: domahidi@control.ee.ethz.ch                                   *\n");
	PRINTTEXT("  ****************************************************************************\n");
	PRINTTEXT("\n\n");
    PRINTTEXT("PROBLEM SUMMARY:\n");
    PRINTTEXT("    Primal variables (n): %d\n", n);
	PRINTTEXT("Equality constraints (p): %d\n", p);
	PRINTTEXT("     Conic variables (m): %d\n", m);
	PRINTTEXT("- - - - - - - - - - - - - - -\n");
    PRINTTEXT("         Size of LP cone: %d\n", l);
    PRINTTEXT("          Number of SOCs: %d\n", ncones);
    for( i=0; i<ncones; i++ ){
        PRINTTEXT("    Size of SOC #%02d: %d\n", i+1, q[i]);
    }
#endif
	
	/* get work data structure */
    mywork = (pwork *)MALLOC(sizeof(pwork));
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for WORK struct\n");
#endif

	/* dimensions */
	mywork->n = n;
	mywork->m = m;
	mywork->p = p;
#if PRINTLEVEL > 2
    PRINTTEXT("Set dimensions\n");
#endif

	/* variables */
    mywork->x = (pfloat *)MALLOC(n*sizeof(pfloat));
    mywork->y = (pfloat *)MALLOC(p*sizeof(pfloat));
    mywork->z = (pfloat *)MALLOC(m*sizeof(pfloat));
    mywork->s = (pfloat *)MALLOC(m*sizeof(pfloat));
	mywork->lambda = (pfloat *)MALLOC(m*sizeof(pfloat));
	mywork->dsaff_by_W = (pfloat *)MALLOC(m*sizeof(pfloat));
    mywork->dsaff = (pfloat *)MALLOC(m*sizeof(pfloat));
    mywork->dzaff = (pfloat *)MALLOC(m*sizeof(pfloat));
    mywork->saff = (pfloat *)MALLOC(m*sizeof(pfloat));
    mywork->zaff = (pfloat *)MALLOC(m*sizeof(pfloat));
	mywork->W_times_dzaff = (pfloat *)MALLOC(m*sizeof(pfloat));
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for variables\n");
#endif

	/* cones */
	mywork->C = (cone *)MALLOC(sizeof(cone));
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for cone struct\n");
#endif

	/* LP cone */
	mywork->C->lpc = (lpcone *)MALLOC(sizeof(lpcone));
	mywork->C->lpc->p = l;
	if( l > 0 ){
		mywork->C->lpc->w = (pfloat *)MALLOC(l*sizeof(pfloat));
		mywork->C->lpc->v = (pfloat *)MALLOC(l*sizeof(pfloat));
		mywork->C->lpc->kkt_idx = (idxint *)MALLOC(l*sizeof(idxint));
#if PRINTLEVEL > 2
        PRINTTEXT("Memory allocated for LP cone\n");
#endif
	} else {
		mywork->C->lpc->w = NULL;
		mywork->C->lpc->v = NULL;
		mywork->C->lpc->kkt_idx = NULL;
#if PRINTLEVEL > 2
        PRINTTEXT("No LP cone present, pointers filled with NULL\n");
#endif
	}


	/* Second-order cones */
	mywork->C->soc = (socone *)MALLOC(ncones*sizeof(socone));
	mywork->C->nsoc = ncones;
    cidx = 0;
    for( i=0; i<ncones; i++ ){
        conesize = (idxint)q[i];
        mywork->C->soc[i].p = conesize;
        mywork->C->soc[i].szidx = l + cidx;
        mywork->C->soc[i].a = 0;
		mywork->C->soc[i].eta = 0;
        mywork->C->soc[i].etasqrt = 0;
		mywork->C->soc[i].omega = 0;
		mywork->C->soc[i].atilde = 0;
        mywork->C->soc[i].q = (pfloat *)MALLOC((conesize-1)*sizeof(pfloat));
		mywork->C->soc[i].qtilde = (pfloat *)MALLOC((conesize-1)*sizeof(pfloat));
		mywork->C->soc[i].vtilde = (pfloat *)MALLOC((conesize-1)*sizeof(pfloat));
		mywork->C->soc[i].kkt_vU = 0;
		mywork->C->soc[i].kkt_qL = 0;
		mywork->C->soc[i].kkt_qtU = (idxint *)MALLOC((conesize-1)*sizeof(idxint));
        mywork->C->soc[i].sidx = (idxint *)MALLOC((conesize)*sizeof(idxint));
		mywork->C->soc[i].skbar = (pfloat *)MALLOC((conesize)*sizeof(pfloat));
		mywork->C->soc[i].zkbar = (pfloat *)MALLOC((conesize)*sizeof(pfloat));
        cidx += conesize;
    }
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for second-order cones\n");
#endif

	/* info struct */
    mywork->info = (stats *)MALLOC(sizeof(stats));   
	mywork->info->tfactor = 0;
	mywork->info->tkktsolve = 0;
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for info struct\n");
#endif

	/* settings */
	mywork->stgs = (settings *)MALLOC(sizeof(settings));
	mywork->stgs->maxit = MAXIT;
	mywork->stgs->gamma = GAMMA;	
	mywork->stgs->delta = DELTA;	
	//mywork->stgs->nitref = NITREF;
	mywork->stgs->abstol = ABSTOL;	
	mywork->stgs->feastol = FEASTOL;
	mywork->stgs->reltol = RELTOL;
#if PRINTLEVEL > 2
    PRINTTEXT("Written settings\n");
#endif
    
	/* Store problem data */
    mywork->A = createSparseMatrix(p, n, Ajc[n], Ajc, Air, Apr);
	mywork->G = createSparseMatrix(m, n, Gjc[n], Gjc, Gir, Gpr);	
#if PROFILING > 1
	mywork->info->ttranspose = 0;
	tic(&tmattranspose);
#endif
	At = transposeSparseMatrix(mywork->A);
#if PROFILING > 1	
	mywork->info->ttranspose += toc(&tmattranspose);
#endif
#if PRINTLEVEL > 2
    PRINTTEXT("Transposed A\n");
#endif
    
    
#if PROFILING > 1	
	tic(&tmattranspose);
#endif
	Gt = transposeSparseMatrix(mywork->G);    	
#if PROFILING > 1	
	mywork->info->ttranspose += toc(&tmattranspose);
#endif
#if PRINTLEVEL > 2
    PRINTTEXT("Transposed G\n");
#endif
    
    mywork->c = c;
    mywork->h = h;
    mywork->b = b;
#if PRINTLEVEL > 2
    PRINTTEXT("Hung pointers for c, h and b into WORK struct\n");
#endif

	/* size of KKT system */
	nK = n + p + m + ncones;
#if PRINTLEVEL > 2
    PRINTTEXT("Size of KKT system: %d\n", nK);
#endif

	/* KKT system stuff (L comes later after symbolic factorization) */
	mywork->KKT = (kkt *)MALLOC(sizeof(kkt));		
	mywork->KKT->D = (pfloat *)MALLOC(nK*sizeof(pfloat));
	mywork->KKT->Parent = (idxint *)MALLOC(nK*sizeof(idxint));
	mywork->KKT->Pinv = (idxint *)MALLOC(nK*sizeof(idxint));
	mywork->KKT->work1 = (pfloat *)MALLOC(nK*sizeof(pfloat));
	mywork->KKT->work2 = (pfloat *)MALLOC(nK*sizeof(pfloat));
    mywork->KKT->work3 = (pfloat *)MALLOC(nK*sizeof(pfloat));
    mywork->KKT->work4 = (pfloat *)MALLOC(nK*sizeof(pfloat));
    mywork->KKT->work5 = (pfloat *)MALLOC(nK*sizeof(pfloat));
    mywork->KKT->work6 = (pfloat *)MALLOC(nK*sizeof(pfloat));
	mywork->KKT->Flag = (idxint *)MALLOC(nK*sizeof(idxint));	
	mywork->KKT->Pattern = (idxint *)MALLOC(nK*sizeof(idxint));
	mywork->KKT->Sign = (idxint *)MALLOC(nK*sizeof(idxint));
	mywork->KKT->Lnz = (idxint *)MALLOC(nK*sizeof(idxint));	
	mywork->KKT->RHS1 = (pfloat *)MALLOC(nK*sizeof(pfloat));
	mywork->KKT->RHS2 = (pfloat *)MALLOC(nK*sizeof(pfloat));
	mywork->KKT->dx1 = (pfloat *)MALLOC(mywork->n*sizeof(pfloat));
	mywork->KKT->dx2 = (pfloat *)MALLOC(mywork->n*sizeof(pfloat));
	mywork->KKT->dy1 = (pfloat *)MALLOC(mywork->p*sizeof(pfloat));
	mywork->KKT->dy2 = (pfloat *)MALLOC(mywork->p*sizeof(pfloat));
	mywork->KKT->dz1 = (pfloat *)MALLOC(mywork->m*sizeof(pfloat));
	mywork->KKT->dz2 = (pfloat *)MALLOC(mywork->m*sizeof(pfloat));
#if PRINTLEVEL > 2
    PRINTTEXT("Created memory for KKT system\n");
#endif

	/* temporary sign storage */
	Sign = (idxint *)MALLOC(nK*sizeof(idxint));
	DiagIdx = (idxint *)MALLOC((n+p)*sizeof(idxint));

	/* set up KKT system */
#if PROFILING > 1
	tic(&tcreatekkt);
#endif	
	KU = createKKT_U(Gt, At, mywork->stgs->delta, mywork->C, Sign, DiagIdx);
#if PROFILING > 1	
	mywork->info->tkktcreate = toc(&tcreatekkt);
#endif	
#if PRINTLEVEL > 2
    PRINTTEXT("Created upper part of KKT matrix K\n");
#endif

	/* allocate memory in KKT system */
	mywork->KKT->PKPt = newSparseMatrix(nK, nK, KU->nnz);
	mywork->KKT->PK = (idxint *)MALLOC(KU->nnz*sizeof(idxint));		

	/* calculate ordering of KKT matrix */
	P = (idxint *)MALLOC(nK*sizeof(idxint));
#if PROFILING > 1
	tic(&tordering);
#endif
	AMD_defaults(Control);	
	amd_result = AMD_order(nK, KU->jc, KU->ir, P, Control, Info);	
#if PROFILING > 1	
	mywork->info->torder = toc(&tordering);
#endif
#if PRINTLEVEL > 2
	if( amd_result == AMD_OK ){
		PRINTTEXT("AMD ordering successfully computed.\n");
		AMD_info(Info);
	} else {
		PRINTTEXT("Problem in AMD ordering, exiting.\n");
        /* AMD_info(Info); */
        return NULL;
	}	
#endif

	
//#ifndef WITH_STATIC_REG
    /* switch off static regularization */
	//for( i=0; i<n+p; i++ ) { KU->pr[DiagIdx[i]] = 0; }
//#if PRINTLEVEL > 2
  //  PRINTTEXT("Static regularization is switched off");
//#endif
//#else
//#if PRINTLEVEL > 2
//    PRINTTEXT("Static regularization is switched on");
//#endif
//#endif
	
	/* calculate inverse permutation and permutation mapping of KKT matrix */	
	pinv(nK, P, mywork->KKT->Pinv);		
	Pinv = mywork->KKT->Pinv;
#if DEBUG > 0
    dumpDenseMatrix_i(P, nK, 1, "P.txt");
    dumpDenseMatrix_i(mywork->KKT->Pinv, nK, 1, "PINV.txt");
#endif
	permuteSparseSymmetricMatrix(KU, mywork->KKT->Pinv, mywork->KKT->PKPt, mywork->KKT->PK);
    

	/* permute sign vector */    
#if PRINTLEVEL > 2
    PRINTTEXT("Sign = [");
    for( i=0; i<nK; i++ ){ PRINTTEXT("%+d ", Sign[i]); }
    PRINTTEXT("];\n");
#endif
	for( i=0; i<nK; i++ ){ mywork->KKT->Sign[i] = Sign[P[i]]; }
    //for( i=0; i<nK; i++ ){ mywork->KKT->Sign[i] = Sign[i]; }
    
	
	/* symbolic factorization */	
	Ljc = (idxint *)MALLOC((nK+1)*sizeof(idxint));	
	LDL_symbolic2(
		mywork->KKT->PKPt->n,   /* A and L are n-by-n, where n >= 0 */
		mywork->KKT->PKPt->jc,  /* input of size n+1, not modified */
		mywork->KKT->PKPt->ir,	 /* input of size nz=Ap[n], not modified */
		Ljc,					 /* output of size n+1, not defined on input */
		mywork->KKT->Parent,	 /* output of size n, not defined on input */
		mywork->KKT->Lnz,		 /* output of size n, not defined on input */
		mywork->KKT->Flag		 /* workspace of size n, not defn. on input or output */
	);
	

	/* assign memory for L */
	lnz = Ljc[nK];
#if PRINTLEVEL > 2
	PRINTTEXT("Nonzeros in L, excluding diagonal: %d\n", lnz) ;
#endif
	Lir = (idxint *)MALLOC(lnz*sizeof(idxint));
	Lpr = (pfloat *)MALLOC(lnz*sizeof(pfloat));
	mywork->KKT->L = createSparseMatrix(nK, nK, lnz, Ljc, Lir, Lpr);

	/* initialize (3,3) block V to I in KKT matrix */	
	prepareKKT4initU(KU, mywork->C);

	/* permute KKT matrix - we work on this one from now on */
	permuteSparseSymmetricMatrix(KU, mywork->KKT->Pinv, mywork->KKT->PKPt, NULL);

	/* set up RHSp for initialization */
	k = 0; j = 0;
	for( i=0; i<n; i++ ){ mywork->KKT->RHS1[Pinv[k++]] = 0; }
	for( i=0; i<p; i++ ){ mywork->KKT->RHS1[Pinv[k++]] = b[i]; }
	for( i=0; i<l; i++ ){ mywork->KKT->RHS1[Pinv[k++]] = h[i]; j++; }
	for( l=0; l<ncones; l++ ){ 
		for( i=0; i < mywork->C->soc[l].p; i++ ){ mywork->KKT->RHS1[Pinv[k++]] = h[j++]; }
		mywork->KKT->RHS1[Pinv[k++]] = 0;
	}	
	
	/* set up RHSd for initialization */
	for( i=0; i<n; i++ ){ mywork->KKT->RHS2[Pinv[i]] = -c[i]; }
	for( i=n; i<nK; i++ ){ mywork->KKT->RHS2[Pinv[i]] = 0; }

	/* get scalings of problem data */
	rx = norm2(c, n); mywork->resx0 = MAX(1, rx);
	ry = norm2(b, p); mywork->resy0 = MAX(1, ry);
	rz = norm2(h, m); mywork->resz0 = MAX(1, rz);

	/* get memory for residuals */
	mywork->rx = (pfloat *)MALLOC(n*sizeof(pfloat));
	mywork->ry = (pfloat *)MALLOC(p*sizeof(pfloat));
	mywork->rz = (pfloat *)MALLOC(m*sizeof(pfloat));
	
	mywork->info->tsetup = toc(&tsetup);
	
	/* clean up */
	FREE(P);
	freeSparseMatrix(At);
	freeSparseMatrix(Gt);
	freeSparseMatrix(KU);
	FREE(Sign);
	FREE(DiagIdx);

    return mywork;
}
