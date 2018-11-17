/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2012-2015 A. Domahidi [domahidi@embotech.com],
 * Automatic Control Lab, ETH Zurich & embotech GmbH, Zurich, Switzerland.
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
#include "equil.h"

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
#define CALLOC mxCalloc
#else
#define MALLOC malloc
#define FREE free
#define CALLOC calloc
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
 * INPUT:      spmat* Gt - pointer to G'
 *             spmat* At - pointer to A'
 *               cone* C - pointer to cone struct
 *
 * OUTPUT:  idxint* Sign - pointer to vector of signs for regularization
 *              spmat* K - pointer to unpermuted upper triangular part of KKT matrix
 *         idxint* AttoK - vector of indices such that K[AtoK[i]] = A[i]
 *         idxint* GttoK - vector of indices such that K[GtoK[i]] = G[i]
 */
void createKKT_U(spmat* Gt, spmat* At, cone* C, idxint** S, spmat** K,
                 idxint* AttoK, idxint* GttoK)
{
	idxint i, j, k, l, r, row_stop, row, cone_strt, ks, conesize;
	idxint n = Gt->m;
	idxint m = Gt->n;
	idxint p = At ? At->n : 0;
	idxint nK, nnzK;
	pfloat *Kpr = NULL;
	idxint *Kjc = NULL, *Kir = NULL;
    idxint *Sign;
#ifdef EXPCONE
	idxint exp_cone_strt;
#endif

	/* Dimension of KKT matrix
     *   =   n (number of variables)
     *     + p (number of equality constraints)
     *     + m (number of inequality constraints)
     *     + 2*C->nsoc (expansion of SOC scalings)
     */
    nK = n + p + m;
#if CONEMODE == 0
    nK += 2*C->nsoc;
#endif

    /* Number of non-zeros in KKT matrix
     *   =   At->nnz (nnz of equality constraint matrix A)
     *     + Gt->nnz (nnz of inequality constraint matrix)
     *     + C->lpc.p (nnz of LP cone)
     *     + 3*[sum(C->soc[i].p)+1] (nnz of expanded soc scalings)
     *     + 6*C->nexc (3x3 hessian of the exponential cone)
     */
	nnzK = (At ? At->nnz : 0) + Gt->nnz + C->lpc->p;
#if STATICREG == 1
    nnzK += n+p;
#endif
	for( i=0; i<C->nsoc; i++ ){
#if CONEMODE == 0
		nnzK += 3*C->soc[i].p+1;
#else
        nnzK += (C->soc[i].p*(C->soc[i].p+1))/2;
#endif
	}

#ifdef EXPCONE
    nnzK += 6*C->nexc;
#endif


#if PRINTLEVEL > 2
    PRINTTEXT("Non-zeros in KKT matrix: %d\n", (int) nnzK);
#endif

	/* Allocate memory for KKT matrix */
	Kpr = (pfloat *)MALLOC(nnzK*sizeof(pfloat));
	Kir = (idxint *)MALLOC(nnzK*sizeof(idxint));
	Kjc = (idxint *)MALLOC((nK+1)*sizeof(idxint));

    /* Allocate memory for sign vector */
    Sign = (idxint *)MALLOC(nK*sizeof(idxint));
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for sign vector\n");
#endif

	/* Set signs for regularization of (1,1) block */
    for( ks=0; ks < n; ks++ ){
        Sign[ks] = +1; /* (1,1) block */
    }
    for( ks=n; ks < n+p; ks++){
        Sign[ks] = -1; /* (2,2) block */
    }
#if CONEMODE == 0
    for( ks=n+p; ks < n+p+C->lpc->p; ks++){
        Sign[ks] = -1; /* (3,3) block: LP cone */
    }
    ks = n+p+C->lpc->p;
    for( l=0; l<C->nsoc; l++){
        for (i=0; i<C->soc[l].p; i++) {
            Sign[ks++] = -1; /* (3,3) block: SOC, D */
        }
        Sign[ks++] = -1;     /* (3,3) block: SOC, v */
        Sign[ks++] = +1;     /* (3,3) block: SOC, u */
    }
#else
    for (ks=n+p; ks < n+p+m; ks++) {
        Sign[ks] = -1;      /* (3,3) block has -1 sign if all dense */
    }
#endif
#ifdef EXPCONE /* In the case CONEMODE != 0 this is marginally inneficient */
    /* Set the signs of the regularization of the expcone block */
    for(ks=nK-3*C->nexc;ks<nK;)
    {
        Sign[ks++] = -1;
    }
#endif


#if DEBUG > 0
    if (ks!=nK) {
        PRINTTEXT("ks = %d, whereas nK = %d - exiting.", (int)ks, (int)nK);
        exit(-1);
    }
#endif

    /* count the number of non-zero entries in K */
    k = 0;

    /* (1,1) block: the first n columns are empty */
#if STATICREG == 0
    for (j=0; j<n; j++) {
        Kjc[j] = 0;
    }
#else
    for (j=0; j<n; j++) {
        Kjc[j] = j;
        Kir[j] = j;
        Kpr[k++] = DELTASTAT;
    }
#endif

    /* Fill upper triangular part of K with values */
    /* (1,2) block: A' */
	i = 0; /* counter for non-zero entries in A or G, respectively */
	for( j=0; j<p; j++ ){
        /* A' */
		row = At->jc[j];
		row_stop = At->jc[j+1];
		if( row <= row_stop ){
			Kjc[n+j] = k;
			while( row++ < row_stop ){
				Kir[k] = At->ir[i];
				Kpr[k] = At->pr[i];
				AttoK[i++] = k++;
			}
		}
#if STATICREG == 1
        Kir[k] = n+j;
        Kpr[k++] = -DELTASTAT;
#endif
    }
	/* (1,3) and (3,3) block: [G'; 0; -Vinit]
     * where
     *
     *   Vinit = blkdiag(I, blkdiag(I,1,-1), ...,  blkdiag(I,1,-1));
     *                        ^ #number of second-order cones ^
     *
     * Note that we have to prepare the (3,3) block accordingly
     * (put zeros for init but store indices that are used in KKT_update
     * of cone module)
     */

	/* LP cone */
	i = 0;
	for( j=0; j < C->lpc->p; j++ ){
        /* copy in G' */
		row = Gt->jc[j];
		row_stop = Gt->jc[j+1];
		if( row <= row_stop ){
			Kjc[n+p+j] = k;
			while( row++ < row_stop ){
				Kir[k] = Gt->ir[i];
				Kpr[k] = Gt->pr[i];
				GttoK[i++] = k++;
			}
		}
        /* -I for LP-cone */
		C->lpc->kkt_idx[j] = k;
		Kir[k] = n+p+j;
		Kpr[k] = -1.0;
        k++;
	}

#if CONEMODE == 0
	/* Second-order cones - copy in G' and set up the scaling matrix
     * which has the following structure:
     *
     *
     *    [ *             0  * ]
     *    [   *           *  * ]
     *    [     *         *  * ]       [ I   v  u  ]      I: identity of size conesize
     *  - [       *       *  * ]   =  -[ u'  1  0  ]      v: vector of size conesize - 1
     *    [         *     *  * ]       [ v'  0' -1 ]      u: vector of size conesize
     *    [           *   *  * ]
     *    [             * *  * ]
     *    [ 0 * * * * * * 1  0 ]
     *    [ * * * * * * * 0 -1 ]
     *
     * NOTE: only the upper triangular part (with the diagonal elements)
     *       is copied in here.
     */
	cone_strt = C->lpc->p;
    for( l=0; l < C->nsoc; l++ ){

        /* size of the cone */
        conesize = C->soc[l].p;

        /* go column-wise about it */
		for( j=0; j < conesize; j++ ){

           row = Gt->jc[cone_strt+j];
           row_stop = Gt->jc[cone_strt+j+1];
           if( row <= row_stop ){
               Kjc[n+p+cone_strt+2*l+j] = k;
               while( row++ < row_stop ){
                   Kir[k] = Gt->ir[i];
                   Kpr[k] = Gt->pr[i];
                   GttoK[i++] = k++;
               }
           }

           /* diagonal D */
           Kir[k] = n+p+cone_strt + 2*l + j;
           Kpr[k] = -1.0;
           C->soc[l].Didx[j] = k;
           k++;
        }

        /* v */
        Kjc[n+p+cone_strt+2*l+conesize] = k;
        for (r=1; r<conesize; r++) {
            Kir[k] = n+p+cone_strt + 2*l + r;
            Kpr[k] = 0;
            k++;
        }
        Kir[k] = n+p+cone_strt + 2*l + conesize;
        Kpr[k] = -1;
        k++;


        /* u */
        Kjc[n+p+cone_strt+2*l+conesize+1] = k;
        for (r=0; r<conesize; r++) {
            Kir[k] = n+p+cone_strt + 2*l + r;
            Kpr[k] = 0;
            k++;
        }
        Kir[k] = n+p+cone_strt + 2*l + conesize + 1;
        Kpr[k] = +1;
        k++;


        /* prepare index for next cone */
		cone_strt += C->soc[l].p;
	}
#else
    /* Second-order cones - copy in G' and set up the scaling matrix
     * which has a dense structure (only upper half part is shown):
     *
     *                     [ a        b*q'     ]
     *  - W^2 = -V = eta^2 [ b*q  I + c*(q*q') ]
     *
     * where    I: identity of size conesize
     *          q: vector of size consize - 1
     *      a,b,c: scalars
     *
     * NOTE: only the upper triangular part (with the diagonal elements)
     *       is copied in here.
     */
	cone_strt = C->lpc->p;
    for( l=0; l < C->nsoc; l++ ){

        /* size of the cone */
        conesize = C->soc[l].p;

        /* go column-wise about it */
		for( j=0; j < conesize; j++ ){

            row = Gt->jc[cone_strt+j];
            row_stop = Gt->jc[cone_strt+j+1];
            if( row <= row_stop ){
                Kjc[n+p+cone_strt+j] = k;
                while( row++ < row_stop ){
                    Kir[k] = Gt->ir[i];
                    Kpr[k] = Gt->pr[i];
                    GttoK[i++] = k++;
                }
            }

            /* first elements - record where this column starts */
            Kir[k] = n+p+cone_strt;
            Kpr[k] = -1.0;
            C->soc[l].colstart[j] = k;
            k++;

            /* the rest of the column */
            for (r=1; r<=j; r++) {
                Kir[k] = n+p+cone_strt+r;
                Kpr[k] = -1.0;
                k++;
            }
        }

        /* prepare index for next cone */
		cone_strt += C->soc[l].p;
	}
#endif
#ifdef EXPCONE
    /* Exponential cones - copy in G' and set up the scaling matrix
     * which has a dense structure (only upper half part is shown):
     *
     *  k is the index of the entry in Kip, Kir
     *  If the dense form of the socp cones is used the exponential cones start
     *  at row n+p+linear_vars+socp_vars.
     *  If the sparse form of the socp cones is used then the exponential cones will
     *  start at row n+p+linear_vars+socp_vars+2*number_socp_cones
     *  In either case cone_strt will have value equal to socp_vars+linear_vars.
     */

    /* At this point cone_strt = sum(q)+l */

    exp_cone_strt = cone_strt;
#if CONEMODE == 0
    exp_cone_strt += 2*C->nsoc;
#endif
     for( l=0; l < C->nexc; l++ ){

        /* go column-wise about it */
		for( j=0; j < 3; j++ ){
            row = Gt->jc[cone_strt+j];
            row_stop = Gt->jc[cone_strt+j+1];
            if( row <= row_stop ){
                Kjc[n+p+exp_cone_strt+j] = k;
                while( row++ < row_stop ){
                    Kir[k] = Gt->ir[i];
                    Kpr[k] = Gt->pr[i];
                    GttoK[i++] = k++;
                }
            }

            /* Save the index of the start of the column for the permutation */
            C->expc[l].colstart[j] = k;
            /* top of the column */
            for (r=0; r<j; r++) {
                Kir[k] = n+p+exp_cone_strt+r;
                Kpr[k] = 0.0;
                k++;
            }
            /*diagonal*/
            Kir[k]   = n+p+exp_cone_strt+j;
            Kpr[k++] = -1.0;
        }

        /* prepare index for next cone */
		cone_strt += 3;
        exp_cone_strt += 3;
	}
#endif


#if PRINTLEVEL > 2
    PRINTTEXT("CREATEKKT: Written %d KKT entries\n", (int)k);
    PRINTTEXT("CREATEKKT: nK=%d and ks=%d\n",(int)nK,(int)ks);
    PRINTTEXT("CREATEKKT: Size of KKT matrix: %d\n", (int)nK);
#endif

	/* return Sign vector and KKT matrix */
    *S = Sign;
	*K = ecoscreateSparseMatrix(nK, nK, nnzK, Kjc, Kir, Kpr);
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

#if defined EQUILIBRATE && EQUILIBRATE > 0
    /* unset equilibration */
    unset_equilibration(w);
#endif

	/* Free KKT related memory      ---            below are the corresponding MALLOCs                */
	FREE(w->KKT->D);                /* mywork->KKT->D = (pfloat *)MALLOC(nK*sizeof(pfloat));          */
	FREE(w->KKT->dx1);              /* mywork->KKT->dx1 = (pfloat *)MALLOC(mywork->n*sizeof(pfloat)); */
	FREE(w->KKT->dx2);              /* mywork->KKT->dx2 = (pfloat *)MALLOC(mywork->n*sizeof(pfloat)); */
	FREE(w->KKT->dy1);              /* mywork->KKT->dy1 = (pfloat *)MALLOC(mywork->p*sizeof(pfloat)); */
	FREE(w->KKT->dy2);              /* mywork->KKT->dy2 = (pfloat *)MALLOC(mywork->p*sizeof(pfloat)); */
	FREE(w->KKT->dz1);              /* mywork->KKT->dz1 = (pfloat *)MALLOC(mywork->m*sizeof(pfloat)); */
	FREE(w->KKT->dz2);              /* mywork->KKT->dz2 = (pfloat *)MALLOC(mywork->m*sizeof(pfloat)); */
	FREE(w->KKT->Flag);             /* mywork->KKT->Flag = (idxint *)MALLOC(nK*sizeof(idxint));       */
	freeSparseMatrix(w->KKT->L);
	FREE(w->KKT->Lnz);              /* mywork->KKT->Lnz = (idxint *)MALLOC(nK*sizeof(idxint));        */
	FREE(w->KKT->Parent);           /* mywork->KKT->Parent = (idxint *)MALLOC(nK*sizeof(idxint));     */
	FREE(w->KKT->Pattern);          /* mywork->KKT->Pattern = (idxint *)MALLOC(nK*sizeof(idxint));    */
	FREE(w->KKT->Sign);             /* mywork->KKT->Sign = (idxint *)MALLOC(nK*sizeof(idxint));       */
	FREE(w->KKT->Pinv);             /* mywork->KKT->Pinv = (idxint *)MALLOC(nK*sizeof(idxint));       */
    FREE(w->KKT->P);
	FREE(w->KKT->PK);               /* mywork->KKT->PK = (idxint *)MALLOC(KU->nnz*sizeof(idxint));    */
	freeSparseMatrix(w->KKT->PKPt); /* mywork->KKT->PKPt = newSparseMatrix(nK, nK, KU->nnz);          */
	FREE(w->KKT->RHS1);             /* mywork->KKT->RHS1 = (pfloat *)MALLOC(nK*sizeof(pfloat));       */
	FREE(w->KKT->RHS2);             /* mywork->KKT->RHS2 = (pfloat *)MALLOC(nK*sizeof(pfloat));       */
	FREE(w->KKT->work1);            /* mywork->KKT->work1 = (pfloat *)MALLOC(nK*sizeof(pfloat));      */
	FREE(w->KKT->work2);            /* mywork->KKT->work2 = (pfloat *)MALLOC(nK*sizeof(pfloat));      */
    FREE(w->KKT->work3);            /* mywork->KKT->work3 = (pfloat *)MALLOC(nK*sizeof(pfloat));      */
    FREE(w->KKT->work4);            /* mywork->KKT->work4 = (pfloat *)MALLOC(nK*sizeof(pfloat));      */
    FREE(w->KKT->work5);            /* mywork->KKT->work5 = (pfloat *)MALLOC(nK*sizeof(pfloat));      */
    FREE(w->KKT->work6);            /* mywork->KKT->work6 = (pfloat *)MALLOC(nK*sizeof(pfloat));      */
	FREE(w->KKT);                   /* mywork->KKT = (kkt *)MALLOC(sizeof(kkt));                      */
	if (w->A) {
		FREE(w->AtoK);
	}
	FREE(w->GtoK);

	/* Free memory for cones */
	if( w->C->lpc->p > 0 ){
		FREE(w->C->lpc->kkt_idx);
		FREE(w->C->lpc->v);
		FREE(w->C->lpc->w);
	}
    /* C->lpc is always allocated, so we free it here. */
    FREE(w->C->lpc);

	for( i=0; i < w->C->nsoc; i++ ){
		FREE(w->C->soc[i].q);
		FREE(w->C->soc[i].skbar);
		FREE(w->C->soc[i].zkbar);
#if CONEMODE == 0
		FREE(w->C->soc[i].Didx);
#endif
#if CONEMODE > 0
        FREE(w->C->soc[i].colstart);
#endif
	}
	if( w->C->nsoc > 0 ){
		FREE(w->C->soc);
	}
#ifdef EXPCONE
    if(w->C->nexc > 0){
        FREE(w->C->expc);
    }
#endif
	FREE(w->C);

	/* free stuff from pwork */
    FREE(w->W_times_dzaff);
	FREE(w->dsaff_by_W);
    FREE(w->dzaff);
    FREE(w->dsaff);
    FREE(w->zaff);
    FREE(w->saff);
	FREE(w->info);
    FREE(w->best_info);
	FREE(w->lambda);
	FREE(w->rx);
	FREE(w->ry);
	FREE(w->rz);
	FREE(w->stgs);
	FREE(w->G);
	if( w->A ) FREE(w->A);
    FREE(w->best_z);
    FREE(w->best_s);
    FREE(w->best_y);
    FREE(w->best_x);
	if( keepvars < 4 ) { FREE(w->z); }
	if( keepvars < 3 ) { FREE(w->s); }
	if( keepvars < 2 ) { FREE(w->y); }
	if( keepvars < 1 ) { FREE(w->x); }
#if defined EQUILIBRATE && EQUILIBRATE > 0
    FREE(w->xequil);
    FREE(w->Aequil);
    FREE(w->Gequil);
#endif
	FREE(w);
}


/*
 * Sets up all data structures needed.
 * Replace by codegen
 */
pwork* ECOS_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones, idxint* q, idxint nexc,
                   pfloat* Gpr, idxint* Gjc, idxint* Gir,
                   pfloat* Apr, idxint* Ajc, idxint* Air,
                   pfloat* c, pfloat* h, pfloat* b)
{
    idxint i, cidx, conesize, lnz, amd_result, nK, *Ljc, *Lir, *P, *Pinv, *Sign;
    pwork* mywork;
	double Control [AMD_CONTROL], Info [AMD_INFO];
	pfloat *Lpr;
	spmat *At, *Gt, *KU;
	idxint *AtoAt, *GtoGt, *AttoK, *GttoK;

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
	PRINTTEXT("  *******************************************************************************\n");
	PRINTTEXT("  * ECOS: Embedded Conic Solver - Sparse Interior Point method for SOCPs        *\n");
	PRINTTEXT("  *                                                                             *\n");
	PRINTTEXT("  * NOTE: The solver is based on L. Vandenberghe's 'The CVXOPT linear and quad- *\n");
	PRINTTEXT("  *       ratic cone program solvers', March 20, 2010. Available online:        *\n");
	PRINTTEXT("  *       [http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf]           *\n");
	PRINTTEXT("  *                                                                             *\n");
	PRINTTEXT("  *       This code uses T.A. Davis' sparse LDL package and AMD code.           *\n");
	PRINTTEXT("  *       [http://www.cise.ufl.edu/research/sparse]                             *\n");
	PRINTTEXT("  *                                                                             *\n");
	PRINTTEXT("  *       Written during a summer visit at Stanford University with S. Boyd.    *\n");
	PRINTTEXT("  *                                                                             *\n");
	PRINTTEXT("  * (C) Alexander Domahidi, ETH Zurich & embotech GmbH, Switzerland, 2012-14.   *\n");
	PRINTTEXT("  *                     Email: domahidi@embotech.com                            *\n");
	PRINTTEXT("  *                                                                             *\n");
    PRINTTEXT("  * Exponential cone solver: Santiago Akle, 2015 tiagoakle@gmail.com            *\n");
	PRINTTEXT("  *******************************************************************************\n");
	PRINTTEXT("\n\n");
    PRINTTEXT("PROBLEM SUMMARY:\n");
    PRINTTEXT("    Primal variables (n): %d\n", (int)n);
	PRINTTEXT("Equality constraints (p): %d\n", (int)p);
	PRINTTEXT("     Conic variables (m): %d\n", (int)m);
	PRINTTEXT("- - - - - - - - - - - - - - -\n");
    PRINTTEXT("         Size of LP cone: %d\n", (int)l);
    PRINTTEXT("          Number of SOCs: %d\n", (int)ncones);
    for( i=0; i<ncones; i++ ){
        PRINTTEXT("    Size of SOC #%02d: %d\n", (int)(i+1), (int)q[i]);
    }
#if defined EXPCONE
    PRINTTEXT("          Number of EXCs: %d\n", (int)nexc);
#endif
    PRINTTEXT("\n");

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
    mywork->D = l + ncones;
#ifdef EXPCONE
    mywork->D = mywork->D + 3*nexc;
#endif
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

    /* best iterates so far */
    mywork->best_x = (pfloat *)MALLOC(n*sizeof(pfloat));
    mywork->best_y = (pfloat *)MALLOC(p*sizeof(pfloat));
    mywork->best_z = (pfloat *)MALLOC(m*sizeof(pfloat));
    mywork->best_s = (pfloat *)MALLOC(m*sizeof(pfloat));
    mywork->best_info = (stats *)MALLOC(sizeof(stats));

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
	mywork->C->soc = (ncones == 0) ? NULL : (socone *)MALLOC(ncones*sizeof(socone));
	mywork->C->nsoc = ncones;
    cidx = 0;
    for( i=0; i<ncones; i++ ){
        conesize = (idxint)q[i];
        mywork->C->soc[i].p = conesize;
        mywork->C->soc[i].a = 0;
		mywork->C->soc[i].eta = 0;
        mywork->C->soc[i].q = (pfloat *)MALLOC((conesize-1)*sizeof(pfloat));
		mywork->C->soc[i].skbar = (pfloat *)MALLOC((conesize)*sizeof(pfloat));
		mywork->C->soc[i].zkbar = (pfloat *)MALLOC((conesize)*sizeof(pfloat));
#if CONEMODE == 0
        mywork->C->soc[i].Didx = (idxint *)MALLOC((conesize)*sizeof(idxint));
#endif
#if CONEMODE > 0
        mywork->C->soc[i].colstart = (idxint *)MALLOC((conesize)*sizeof(idxint));
#endif
        cidx += conesize;
    }
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for second-order cones\n");
#endif

#ifdef EXPCONE
    /* Exponential cones */
    mywork->C->nexc  = nexc;
    mywork->C->expc  = (nexc == 0) ? NULL : (expcone *)MALLOC(nexc*sizeof(expcone));
    mywork->C->fexv  = cidx+l;
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for exponential cones\n");
#endif
#endif

    /* The number of conic variables has to equal l+sum(q)+3*nexp else terminate*/
#ifdef EXPCONE
    if(cidx+l+3*nexc!=m)
    {
#if PRINTLEVEL > 2
        PRINTTEXT("Number of conic variables does not match l+sum(q)+3*nexc\n");
#endif
        return NULL;
    }
#else
    if(cidx+l!=m)
    {
#if PRINTLEVEL > 2
        PRINTTEXT("Number of conic variables does not match l+sum(q)\n");
#endif
        return NULL;
    }
#endif

	/* info struct */
    mywork->info = (stats *)MALLOC(sizeof(stats));
#if PROFILING > 1
	mywork->info->tfactor = 0;
	mywork->info->tkktsolve = 0;
    mywork->info->tfactor_t1 = 0;
    mywork->info->tfactor_t2 = 0;
#endif
#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for info struct\n");
#endif

#if defined EQUILIBRATE && EQUILIBRATE > 0
    /* equilibration vector */
    mywork->xequil = (pfloat *)MALLOC(n*sizeof(pfloat));
    mywork->Aequil = (pfloat *)MALLOC(p*sizeof(pfloat));
    mywork->Gequil = (pfloat *)MALLOC(m*sizeof(pfloat));

#if PRINTLEVEL > 2
    PRINTTEXT("Memory allocated for equilibration vectors\n");
#endif
#endif

	/* settings */
	mywork->stgs = (settings *)MALLOC(sizeof(settings));
	mywork->stgs->maxit = MAXIT;
	mywork->stgs->gamma = GAMMA;
	mywork->stgs->delta = DELTA;
    mywork->stgs->eps = EPS;
	mywork->stgs->nitref = NITREF;
	mywork->stgs->abstol = ABSTOL;
	mywork->stgs->feastol = FEASTOL;
	mywork->stgs->reltol = RELTOL;
    mywork->stgs->abstol_inacc = ATOL_INACC;
	mywork->stgs->feastol_inacc = FTOL_INACC;
	mywork->stgs->reltol_inacc = RTOL_INACC;
    mywork->stgs->verbose = VERBOSE;
    #ifdef EXPCONE
   	mywork->stgs->max_bk_iter = MAX_BK;
    mywork->stgs->bk_scale    = BK_SCALE;
    mywork->stgs->centrality  = CENTRALITY;
#endif


#if PRINTLEVEL > 2
    PRINTTEXT("Written settings\n");
#endif

    mywork->c = c;
    mywork->h = h;
    mywork->b = b;
#if PRINTLEVEL > 2
    PRINTTEXT("Hung pointers for c, h and b into WORK struct\n");
#endif

    /* Store problem data */
  if(Apr && Ajc && Air) {
    mywork->A = ecoscreateSparseMatrix(p, n, Ajc[n], Ajc, Air, Apr);
  } else {
    mywork->A = NULL;
  }
  if (Gpr && Gjc && Gir) {
	  mywork->G = ecoscreateSparseMatrix(m, n, Gjc[n], Gjc, Gir, Gpr);
  } else {
    /* create an empty sparse matrix */
	mywork->G = ecoscreateSparseMatrix(m, n, 0, Gjc, Gir, Gpr);
  }

#if defined EQUILIBRATE && EQUILIBRATE > 0
    set_equilibration(mywork);
    #if PRINTLEVEL > 2
        PRINTTEXT("Done equilibrating\n");
    #endif
#endif

#if PROFILING > 1
	mywork->info->ttranspose = 0;
	tic(&tmattranspose);
#endif
  if(mywork->A) {
      AtoAt = MALLOC(mywork->A->nnz*sizeof(idxint));
	  At = transposeSparseMatrix(mywork->A, AtoAt);
  }
  else {
    At = NULL;
    AtoAt = NULL;
  }
#if PROFILING > 1
	mywork->info->ttranspose += toc(&tmattranspose);
#endif
#if PRINTLEVEL > 2
    PRINTTEXT("Transposed A\n");
#endif


#if PROFILING > 1
	tic(&tmattranspose);
#endif
	GtoGt = MALLOC(mywork->G->nnz*sizeof(idxint));
	Gt = transposeSparseMatrix(mywork->G, GtoGt);
#if PROFILING > 1
	mywork->info->ttranspose += toc(&tmattranspose);
#endif
#if PRINTLEVEL > 2
    PRINTTEXT("Transposed G\n");
#endif

    /* set up KKT system */
#if PROFILING > 1
	tic(&tcreatekkt);
#endif
	if (mywork->A)
		AttoK = MALLOC(mywork->A->nnz*sizeof(idxint));
	else
		AttoK = NULL;
	GttoK = MALLOC(mywork->G->nnz*sizeof(idxint));
	createKKT_U(Gt, At, mywork->C, &Sign, &KU, AttoK, GttoK);
#if PROFILING > 1
	mywork->info->tkktcreate = toc(&tcreatekkt);
#endif
#if PRINTLEVEL > 2
    PRINTTEXT("Created upper part of KKT matrix K\n");
#endif

    /* Save a mapping from data in A and G to the KKT matrix */
    if (mywork->A){
        mywork->AtoK = MALLOC(mywork->A->nnz*sizeof(idxint));
        for(i=0; i<mywork->A->nnz; i++){ mywork->AtoK[i] = AttoK[AtoAt[i]]; }
    }
    else
        mywork->AtoK = NULL;
    mywork->GtoK = MALLOC(mywork->G->nnz*sizeof(idxint));
    for(i=0; i<mywork->G->nnz; i++){ mywork->GtoK[i] = GttoK[GtoGt[i]]; }


	/*
     * Set up KKT system related data
     * (L comes later after symbolic factorization)
     */
    nK = KU->n;

#if DEBUG > 0
    dumpSparseMatrix(KU, "KU0.txt");
#endif
#if PRINTLEVEL > 2
    PRINTTEXT("Dimension of KKT matrix: %d\n", (int)nK);
    PRINTTEXT("Non-zeros in KKT matrix: %d\n", (int)KU->nnz);
#endif



    /* allocate memory in KKT system */
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
	mywork->KKT->Lnz = (idxint *)MALLOC(nK*sizeof(idxint));
	mywork->KKT->RHS1 = (pfloat *)MALLOC(nK*sizeof(pfloat));
	mywork->KKT->RHS2 = (pfloat *)MALLOC(nK*sizeof(pfloat));
	mywork->KKT->dx1 = (pfloat *)MALLOC(mywork->n*sizeof(pfloat));
	mywork->KKT->dx2 = (pfloat *)MALLOC(mywork->n*sizeof(pfloat));
	mywork->KKT->dy1 = (pfloat *)MALLOC(mywork->p*sizeof(pfloat));
	mywork->KKT->dy2 = (pfloat *)MALLOC(mywork->p*sizeof(pfloat));
	mywork->KKT->dz1 = (pfloat *)MALLOC(mywork->m*sizeof(pfloat));
	mywork->KKT->dz2 = (pfloat *)MALLOC(mywork->m*sizeof(pfloat));
    mywork->KKT->Sign = (idxint *)MALLOC(nK*sizeof(idxint));
    mywork->KKT->PKPt = newSparseMatrix(nK, nK, KU->nnz);
	mywork->KKT->PK = (idxint *)MALLOC(KU->nnz*sizeof(idxint));

#if PRINTLEVEL > 2
    PRINTTEXT("Created memory for KKT-related data\n");
#endif



    /* calculate ordering of KKT matrix using AMD */
	P = (idxint *)MALLOC(nK*sizeof(idxint));
#if PROFILING > 1
	tic(&tordering);
#endif
	AMD_defaults(Control);
	amd_result = AMD_order(nK, KU->jc, KU->ir, P, Control, Info);
#if PROFILING > 1
	mywork->info->torder = toc(&tordering);
#endif

	if( amd_result == AMD_OK ){
#if PRINTLEVEL > 2
		PRINTTEXT("AMD ordering successfully computed.\n");
		AMD_info(Info);
#endif
	} else {
#if PRINTLEVEL > 2
		PRINTTEXT("Problem in AMD ordering, exiting.\n");
        AMD_info(Info);
#endif
        return NULL;
	}

	/* calculate inverse permutation and permutation mapping of KKT matrix */
	pinv(nK, P, mywork->KKT->Pinv);
	Pinv = mywork->KKT->Pinv;
#if DEBUG > 0
    dumpDenseMatrix_i(P, nK, 1, "P.txt");
    dumpDenseMatrix_i(mywork->KKT->Pinv, nK, 1, "PINV.txt");
#endif
	permuteSparseSymmetricMatrix(KU, mywork->KKT->Pinv, mywork->KKT->PKPt, mywork->KKT->PK);

	/* permute sign vector */
    for( i=0; i<nK; i++ ){ mywork->KKT->Sign[Pinv[i]] = Sign[i]; }
#if PRINTLEVEL > 3
    PRINTTEXT("P = [");
    for( i=0; i<nK; i++ ){ PRINTTEXT("%d ", (int)P[i]); }
    PRINTTEXT("];\n");
    PRINTTEXT("Pinv = [");
    for( i=0; i<nK; i++ ){ PRINTTEXT("%d ", (int)Pinv[i]); }
    PRINTTEXT("];\n");
    PRINTTEXT("Sign = [");
    for( i=0; i<nK; i++ ){ PRINTTEXT("%+d ", (int)Sign[i]); }
    PRINTTEXT("];\n");
    PRINTTEXT("SignP = [");
    for( i=0; i<nK; i++ ){ PRINTTEXT("%+d ", (int)mywork->KKT->Sign[i]); }
    PRINTTEXT("];\n");
#endif



	/* symbolic factorization */
	Ljc = (idxint *)MALLOC((nK+1)*sizeof(idxint));
#if PRINTLEVEL > 2
    PRINTTEXT("Allocated memory for cholesky factor L\n");
#endif
	LDL_symbolic2(
		mywork->KKT->PKPt->n,    /* A and L are n-by-n, where n >= 0 */
		mywork->KKT->PKPt->jc,   /* input of size n+1, not modified */
		mywork->KKT->PKPt->ir,	 /* input of size nz=Ap[n], not modified */
		Ljc,					 /* output of size n+1, not defined on input */
		mywork->KKT->Parent,	 /* output of size n, not defined on input */
		mywork->KKT->Lnz,		 /* output of size n, not defined on input */
		mywork->KKT->Flag		 /* workspace of size n, not defn. on input or output */
	);


	/* assign memory for L */
	lnz = Ljc[nK];
#if PRINTLEVEL > 2
	PRINTTEXT("Nonzeros in L, excluding diagonal: %d\n", (int)lnz) ;
#endif
	Lir = (idxint *)MALLOC(lnz*sizeof(idxint));
	Lpr = (pfloat *)MALLOC(lnz*sizeof(pfloat));
	mywork->KKT->L = ecoscreateSparseMatrix(nK, nK, lnz, Ljc, Lir, Lpr);
#if PRINTLEVEL > 2
	PRINTTEXT("Created Cholesky factor of K in KKT struct\n");
#endif


	/* permute KKT matrix - we work on this one from now on */
	permuteSparseSymmetricMatrix(KU, mywork->KKT->Pinv, mywork->KKT->PKPt, NULL);
#if DEBUG > 0
    dumpSparseMatrix(mywork->KKT->PKPt, "PKPt.txt");
#endif

	/* get memory for residuals */
	mywork->rx = (n == 0) ? NULL : (pfloat *)MALLOC(n*sizeof(pfloat));
	mywork->ry = (p == 0) ? NULL : (pfloat *)MALLOC(p*sizeof(pfloat));
	mywork->rz = (m == 0) ? NULL : (pfloat *)MALLOC(m*sizeof(pfloat));

    /* clean up */
    mywork->KKT->P = P;
	FREE(Sign);
    if(At) {
        freeSparseMatrix(At);
        FREE(AtoAt);
        FREE(AttoK);
    }
	freeSparseMatrix(Gt);
	freeSparseMatrix(KU);
    FREE(GtoGt);
    FREE(GttoK);

#if PROFILING > 0
	mywork->info->tsetup = toc(&tsetup);
#endif

    return mywork;
}
