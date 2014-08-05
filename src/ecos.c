/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2012-14 Alexander Domahidi [domahidi@control.ee.ethz.ch],
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

/* Main solver module */


/* ECOS HEADER FILE ---------------------------------------------------- */
#include "ecos.h"

#include "splamm.h"

/* NEEDED FOR SQRT ----------------------------------------------------- */
#include <math.h>


/* Some internal defines */
#define ECOS_NOT_CONVERGED_YET (-87)  /* indicates no convergence yet    */

/**
 * Private static const char * for version numbering.
 * All versions point to this string.
 */
static const char* thisVersion = ECOS_VERSION;


/**
 * Version: returns the current version number
 * Returns the version number in the format X.Y.Z
 */
const char* ECOS_ver(void)
{
    return thisVersion;
}


/* Compares stats of two iterates with each other.
 * Returns 1 if infoA is better than infoB, zero otherwise.
 */
idxint compareStatistics(stats* infoA, stats* infoB)
{
   
    if ( infoA->pinfres != NAN && infoA->kapovert > 1){
        if( infoB->pinfres != NAN ) {
            /* A->pinfres != NAN, B->pinfres!=NAN */
            if ( ( infoA->gap > 0 && infoB->gap > 0 && infoA->gap < infoB->gap ) &&
                ( infoA->pinfres > 0 && infoA->pinfres < infoB->pres ) &&
                ( infoA->mu > 0 && infoA->mu < infoB->mu ) ){
                /* PRINTTEXT("BRANCH 1 "); */
                return 1;
            } else {
                /* PRINTTEXT("BRANCH 1 not OK"); */
                return 0;
            }
        } else {
            /* A->pinfres != NAN, B->pinfres==NAN */
            if ( ( infoA->gap > 0 && infoB->gap > 0 && infoA->gap < infoB->gap ) &&
               ( infoA->mu > 0 && infoA->mu < infoB->mu ) ){
                /* PRINTTEXT("BRANCH 2 "); */
                return 1;
            } else {
                /* PRINTTEXT("BRANCH 2 not OK"); */
                return 0;
            }
        }
    } else {
            /* A->pinfres == NAN or pinfres too large */
        if ( ( infoA->gap > 0 && infoB->gap > 0 && infoA->gap < infoB->gap ) &&
            ( infoA->pres > 0 && infoA->pres < infoB->pres ) &&
            ( infoA->dres > 0 && infoA->dres < infoB->dres ) &&
            ( infoA->kapovert > 0 && infoA->kapovert < infoB->kapovert) &&
            ( infoA->mu > 0 && infoA->mu < infoB->mu ) ){
            /* PRINTTEXT("BRANCH 3 OK "); */
            return 1;
        } else {
            /* PRINTTEXT("BRANCH 3 not OK"); */
            return 0;
        }
    }
}


/* Copy variables from current to best iterate */
void saveIterateAsBest(pwork* w)
{
    idxint i;
    for (i=0; i<w->n; i++) { w->best_x[i] = w->x[i]; }
    for (i=0; i<w->p; i++) { w->best_y[i] = w->y[i]; }
    for (i=0; i<w->m; i++) { w->best_z[i] = w->z[i]; }
    for (i=0; i<w->m; i++) { w->best_s[i] = w->s[i]; }
    w->best_kap = w->kap;
    w->best_tau = w->tau;
    w->best_cx = w->cx;
    w->best_by = w->by;
    w->best_hz = w->hz;
    w->best_info->pcost = w->info->pcost;
	w->best_info->dcost = w->info->dcost;
    w->best_info->pres = w->info->pres;
    w->best_info->dres = w->info->dres;
	w->best_info->pinfres = w->info->pinfres;
    w->best_info->dinfres = w->info->dinfres;
    w->best_info->gap = w->info->gap;
    w->best_info->relgap = w->info->relgap;
    w->best_info->mu = w->info->mu;
	w->best_info->kapovert = w->info->kapovert;
    w->best_info->iter = w->info->iter;
}


/* Copy variables from best iterate to current
 * Comapred to saveIterateAsBest, only the iteration number
 * is not restored.
 */
void restoreBestIterate(pwork* w)
{
    idxint i;
    for (i=0; i<w->n; i++) { w->x[i] = w->best_x[i]; }
    for (i=0; i<w->p; i++) { w->y[i] = w->best_y[i]; }
    for (i=0; i<w->m; i++) { w->z[i] = w->best_z[i]; }
    for (i=0; i<w->m; i++) { w->s[i] = w->best_s[i]; }
    w->kap = w->best_kap;
    w->tau = w->best_tau;
    w->cx = w->best_cx;
    w->by = w->best_by;
    w->hz = w->best_hz;
    w->info->pcost = w->best_info->pcost;
	w->info->dcost = w->best_info->dcost;
    w->info->pres = w->best_info->pres;
    w->info->dres = w->best_info->dres;
	w->info->pinfres = w->best_info->pinfres;
    w->info->dinfres = w->best_info->dinfres;
    w->info->gap = w->best_info->gap;
    w->info->relgap = w->best_info->relgap;
    w->info->mu = w->best_info->mu;
	w->info->kapovert = w->best_info->kapovert;
}


/* 
 * This function is reponsible for checking the exit/convergence conditions of ECOS.
 * If one of the exit conditions is met, ECOS displays an exit message and returns
 * the corresponding exit code. The calling function must then make sure that ECOS
 * is indeed correctly exited, so a call to this function should always be followed
 * by a break statement.
 *
 *    If mode == 0, normal precisions are checked.
 *
 *    If mode != 0, reduced precisions are checked, and the exit display is augmented
 *                  by "Close to". The exitcodes returned are increased by the value
 *                  of mode.
 * 
 * The primal and dual infeasibility flags w->info->pinf and w->info->dinf are raised
 * according to the outcome of the test.
 *
 * If none of the exit tests are met, the function returns ECOS_NOT_CONVERGED_YET.
 * This should not be an exitflag that is ever returned to the outside world.
 */
idxint checkExitConditions(pwork* w, idxint mode)
{
    pfloat feastol;
    pfloat abstol;
    pfloat reltol;
    
    /* set accuracy against which to check */
    if( mode == 0) {
        /* check convergence against normal precisions */
        feastol = w->stgs->feastol;
        abstol = w->stgs->abstol;
        reltol = w->stgs->reltol;
    } else {
        /* check convergence against reduced precisions */
        feastol = w->stgs->feastol_inacc;
        abstol = w->stgs->abstol_inacc;
        reltol = w->stgs->reltol_inacc;
    }
    
    /* Optimal? */
    if( ( -w->cx > 0 || -w->by - w->hz >= -abstol ) &&
        ( w->info->pres < feastol && w->info->dres < feastol ) &&
        ( w->info->gap < abstol || w->info->relgap < reltol  ) ){
#if PRINTLEVEL > 0
        if( w->stgs->verbose ) {
            if( mode == 0) {
                PRINTTEXT("\nOPTIMAL (within feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", MAX(w->info->dres, w->info->pres), w->info->relgap, w->info->gap);
            } else {
                PRINTTEXT("\nClose to OPTIMAL (within feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", MAX(w->info->dres, w->info->pres), w->info->relgap, w->info->gap);
            }
        }
#endif
        w->info->pinf = 0;
        w->info->dinf = 0;
        return ECOS_OPTIMAL + mode;
    }
        
    /* Dual infeasible? */
    else if( (w->info->dinfres != NAN) && (w->info->dinfres < feastol) ){
#if PRINTLEVEL > 0
        if( w->stgs->verbose ) {
            if( mode == 0) {
                PRINTTEXT("\nUNBOUNDED (within feastol=%3.1e).", w->info->dinfres );
            } else {
                PRINTTEXT("\nClose to UNBOUNDED (within feastol=%3.1e).", w->info->dinfres );
            }
        }
#endif
        w->info->pinf = 0;
        w->info->dinf = 1;
        return ECOS_DINF + mode;
    }
    
    /* Primal infeasible? */
    else if( (w->info->pinfres != NAN && w->info->pinfres < feastol) ||
            ( w->tau < w->stgs->feastol && w->kap < w->stgs->feastol && w->info->pinfres < w->stgs->feastol) ){
#if PRINTLEVEL > 0
        if( w->stgs->verbose ) {
            if( mode == 0) {
                PRINTTEXT("\nPRIMAL INFEASIBLE (within feastol=%3.1e).", w->info->pinfres );
            } else {
                PRINTTEXT("\nClose to PRIMAL INFEASIBLE (within feastol=%3.1e).", w->info->pinfres );
            }
        }
#endif
        w->info->pinf = 1;
        w->info->dinf = 0;
        return ECOS_PINF + mode;
    }

    
    /* Indicate if none of the above criteria are met */
    else {
        return ECOS_NOT_CONVERGED_YET;
    }
}


/*
 * Initializes the solver.
 */
idxint init(pwork* w)
{
	idxint i, j, k, l, KKT_FACTOR_RETURN_CODE;
	idxint* Pinv = w->KKT->Pinv;
    pfloat rx, ry, rz;

#if PROFILING > 1
	timer tfactor, tkktsolve;
#endif

    /* set regularization parameter */
	w->KKT->delta = w->stgs->delta;
    
    /* Initialize KKT matrix */
    kkt_init(w->KKT->PKPt, w->KKT->PK, w->C);

#if DEBUG > 0
    dumpSparseMatrix(w->KKT->PKPt, "PKPt0.txt");
#endif

    
    /* initialize RHS1 */
	k = 0; j = 0;
	for( i=0; i<w->n; i++ ){ w->KKT->RHS1[w->KKT->Pinv[k++]] = 0; }
	for( i=0; i<w->p; i++ ){ w->KKT->RHS1[w->KKT->Pinv[k++]] = w->b[i]; }
	for( i=0; i<w->C->lpc->p; i++ ){ w->KKT->RHS1[w->KKT->Pinv[k++]] = w->h[i]; j++; }
	for( l=0; l<w->C->nsoc; l++ ){
		for( i=0; i < w->C->soc[l].p; i++ ){ w->KKT->RHS1[w->KKT->Pinv[k++]] = w->h[j++]; }
#if CONEMODE == 0
		w->KKT->RHS1[w->KKT->Pinv[k++]] = 0;
        w->KKT->RHS1[w->KKT->Pinv[k++]] = 0;
#endif
	}
#if PRINTLEVEL > 2
    PRINTTEXT("Written %d entries of RHS1\n", (int)k);
#endif
	
	/* initialize RHS2 */
	for( i=0; i<w->n; i++ ){ w->KKT->RHS2[w->KKT->Pinv[i]] = -w->c[i]; }
	for( i=w->n; i<w->KKT->PKPt->n; i++ ){ w->KKT->RHS2[w->KKT->Pinv[i]] = 0; }
    
	/* get scalings of problem data */
	rx = norm2(w->c, w->n); w->resx0 = MAX(1, rx);
	ry = norm2(w->b, w->p); w->resy0 = MAX(1, ry);
	rz = norm2(w->h, w->m); w->resz0 = MAX(1, rz);


	/* Factor KKT matrix - this is needed in all 3 linear system solves */
#if PROFILING > 1
	tic(&tfactor);
    KKT_FACTOR_RETURN_CODE = kkt_factor(w->KKT, w->stgs->eps, w->stgs->delta, &w->info->tfactor_t1, &w->info->tfactor_t2);
	w->info->tfactor += toc(&tfactor);
#else
    KKT_FACTOR_RETURN_CODE = kkt_factor(w->KKT, w->stgs->eps, w->stgs->delta);
#endif
    
    /* check if factorization was successful, exit otherwise */
	if(  KKT_FACTOR_RETURN_CODE != KKT_OK ){
#if PRINTLEVEL > 0
    if( w->stgs->verbose ) PRINTTEXT("\nProblem in factoring KKT system, aborting.");
#endif
        return ECOS_FATAL;
    }

	
	/* 
	 * PRIMAL VARIABLES:
	 *  - solve xhat = arg min ||Gx-h||_2^2  such that Ax = b
	 *  - r = h - G*xhat
	 * These two equations are solved by
	 *
	 * [ 0   A'  G' ] [ xhat ]     [ 0 ]
     * [ A   0   0  ] [  y   ]  =  [ b ]
     * [ G   0  -I  ] [ -r   ]     [ h ]
     *
     * and then take shat = r if alphap < 0, zbar + (1+alphap)*e otherwise
     * where alphap = inf{ alpha | sbar + alpha*e >= 0 }
	 */	
	
	/* Solve for RHS [0; b; h] */
#if PROFILING > 1
	tic(&tkktsolve);
#endif
	w->info->nitref1 = kkt_solve(w->KKT, w->A, w->G, w->KKT->RHS1, w->KKT->dx1, w->KKT->dy1, w->KKT->dz1, w->n, w->p, w->m, w->C, 1, w->stgs->nitref);
#if PROFILING > 1
	w->info->tkktsolve += toc(&tkktsolve);
#endif
    
#if DEBUG > 0
#if PRINTLEVEL > 3
    printDenseMatrix(w->KKT->dx1, w->n, 1, "dx1_init");
    printDenseMatrix(w->KKT->dy1, w->p, 1, "dy1_init");
    printDenseMatrix(w->KKT->dz1, w->m, 1, "dz1_init");
#endif
#endif

	/* Copy out initial value of x */
	for( i=0; i<w->n; i++ ){ w->x[i] = w->KKT->dx1[i]; }

	/* Copy out -r into temporary variable */
	for( i=0; i<w->m; i++ ){ w->KKT->work1[i] = -w->KKT->dz1[i]; }

	/* Bring variable to cone */
	bring2cone(w->C, w->KKT->work1, w->s );
	
	/*
	 * dual variables
	 * solve (yhat,zbar) = arg min ||z||_2^2 such that G'*z + A'*y + c = 0
	 *
	 * we can solve this by
	 *
	 * [ 0   A'  G' ] [  x   ]     [ -c ]
	 * [ A   0   0  ] [ yhat ]  =  [  0 ]
	 * [ G   0  -I  ] [ zbar ]     [  0 ]
	 *
	 * and then take zhat = zbar if alphad < 0, zbar + (1+alphad)*e otherwise
	 * where alphad = inf{ alpha | zbar + alpha*e >= 0 }
	 */
	
	/* Solve for RHS [-c; 0; 0] */
#if PROFILING > 1
	tic(&tkktsolve);
#endif
	w->info->nitref2 = kkt_solve(w->KKT, w->A, w->G, w->KKT->RHS2, w->KKT->dx2, w->KKT->dy2, w->KKT->dz2, w->n, w->p, w->m, w->C, 1, w->stgs->nitref);
#if PROFILING > 1
	w->info->tkktsolve += toc(&tkktsolve);
#endif
    
#if DEBUG > 0
#if PRINTLEVEL > 3
    printDenseMatrix(w->KKT->dx2, w->n, 1, "dx2_init");
    printDenseMatrix(w->KKT->dy2, w->p, 1, "dy2_init");
    printDenseMatrix(w->KKT->dz2, w->m, 1, "dz2_init");
#endif
#endif
    
    /* Copy out initial value of y */
	for( i=0; i<w->p; i++ ){ w->y[i] = w->KKT->dy2[i]; }
	
	/* Bring variable to cone */
	bring2cone(w->C, w->KKT->dz2, w->z );
	
	/* Prepare RHS1 - before this line RHS1 = [0; b; h], after it holds [-c; b; h] */
	for( i=0; i<w->n; i++){ w->KKT->RHS1[Pinv[i]] = -w->c[i]; }

	/*
	 * other variables
	 */
	w->kap = 1.0;
	w->tau = 1.0;

	w->info->step = 0;
	w->info->step_aff = 0;
	w->info->dinf = 0;
	w->info->pinf = 0;

	return 0;
}



/*
 * Computes residuals.
 *
 * hrx = -A'*y - G'*z;  rx = hrx - c.*tau;  hresx = norm(rx,2);
 * hry = A*x;           ry = hry - b.*tau;  hresy = norm(ry,2);
 * hrz = s + G*x;       rz = hrz - h.*tau;  hresz = norm(rz,2);
 * rt = kappa + c'*x + b'*y + h'*z;
 */ 
void computeResiduals(pwork *w)
{
	/* rx = -A'*y - G'*z - c.*tau */
	if( w->p > 0 ) {
        sparseMtVm(w->A, w->y, w->rx, 1, 0);
        sparseMtVm(w->G, w->z, w->rx, 0, 0);
    } else {
        sparseMtVm(w->G, w->z, w->rx, 1, 0);
    }
	w->hresx = norm2(w->rx, w->n);
	vsubscale(w->n, w->tau, w->c, w->rx);    
		
	/* ry = A*x - b.*tau */
	if( w->p > 0 ){
        sparseMV(w->A, w->x, w->ry, 1, 1);
        w->hresy = norm2(w->ry, w->p);
        vsubscale(w->p, w->tau, w->b, w->ry);
    } else {
        w->hresy = 0;
        w->ry = NULL;
	}
    
	/* rz = s + G*x - h.*tau */
	sparseMV(w->G, w->x, w->rz, 1, 1);
	vadd(w->m, w->s, w->rz);
	w->hresz = norm2(w->rz, w->m);
	vsubscale(w->m, w->tau, w->h, w->rz);

	/* rt = kappa + c'*x + b'*y + h'*z; */
	w->cx = eddot(w->n, w->c, w->x);
	w->by = w->p > 0 ? eddot(w->p, w->b, w->y) : 0.0;
	w->hz = eddot(w->m, w->h, w->z);
	w->rt = w->kap + w->cx + w->by + w->hz;    
}



/* 
 * Updates statistics.
 */
void updateStatistics(pwork* w)
{
	pfloat nry, nrz;
	
	stats* info = w->info;
	
	/* mu = (s'*z + kap*tau) / (D+1) where s'*z is the duality gap */
	info->gap = eddot(w->m, w->s, w->z);
	info->mu = (info->gap + w->kap*w->tau) / (w->D + 1);

	info->kapovert = w->kap / w->tau;
	info->pcost = w->cx / w->tau;
	info->dcost = -(w->hz + w->by) / w->tau;

	/* relative duality gap */
	if( info->pcost < 0 ){ info->relgap = info->gap / (-info->pcost); }
	else if( info->dcost > 0 ){ info->relgap = info->gap / info->dcost; }
	else info->relgap = NAN;

	/* residuals */
    nry = w->p > 0 ? norm2(w->ry, w->p)/w->resy0 : 0.0;
    nrz = norm2(w->rz, w->m)/w->resz0;
	info->pres = MAX(nry, nrz) / w->tau;
	info->dres = norm2(w->rx, w->n)/w->resx0 / w->tau;
    
	/* infeasibility measures
     *
	 * CVXOPT uses the following:
     * info->pinfres = w->hz + w->by < 0 ? w->hresx / w->resx0 / (-w->hz - w->by) : NAN;
     * info->dinfres = w->cx < 0 ? MAX(w->hresy/w->resy0, w->hresz/w->resz0) / (-w->cx) : NAN;
     */
    info->pinfres = w->hz + w->by < 0 ? w->hresx/w->resx0 : NAN;
    info->dinfres = w->cx < 0 ? MAX(w->hresy/w->resy0, w->hresz/w->resz0) : NAN;
	
    
#if PRINTLEVEL > 2
    PRINTTEXT("TAU=%6.4e  KAP=%6.4e  PINFRES=%6.4e  DINFRES=%6.4e\n",w->tau,w->kap,info->pinfres, info->dinfres );
#endif
    
}



#if PRINTLEVEL > 0
void printProgress(stats* info)
{
	if( info->iter == 0 )
	{
		/* print header at very first iteration */		
#if PRINTLEVEL == 2
		PRINTTEXT("\nECOS %s - (c) A. Domahidi, ETH Zurich & embotech 2012-14. Support: ecos@embotech.com\n\n", ECOS_VERSION);
#endif
#if defined _WIN32 || defined _WIN64		
		PRINTTEXT("It     pcost       dcost      gap   pres   dres    k/t    mu     step    IR\n");
		PRINTTEXT("%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e   N/A    %d %d -\n",(int)info->iter, info->pcost, info->dcost, info->gap, info->pres, info->dres, info->kapovert, info->mu, (int)info->nitref1, (int)info->nitref2);
#else
		PRINTTEXT("It     pcost         dcost      gap     pres    dres     k/t     mu      step     IR\n");
		PRINTTEXT("%2d  %c%+5.3e  %c%+5.3e  %c%+2.0e  %c%2.0e  %c%2.0e  %c%2.0e  %c%2.0e    N/A     %d %d -\n",(int)info->iter, 32, info->pcost, 32, info->dcost, 32, info->gap, 32, info->pres, 32, info->dres, 32, info->kapovert, 32, info->mu, (int)info->nitref1, (int)info->nitref2);
#endif	
	}  else {
#if defined _WIN32 || defined _WIN64
		PRINTTEXT("%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e  %6.4f  %d %d %d\n",(int)info->iter, info->pcost, info->dcost, info->gap, info->pres, info->dres, info->kapovert, info->mu, info->step, (int)info->nitref1, (int)info->nitref2, (int)info->nitref3);
#else
		PRINTTEXT("%2d  %c%+5.3e%c  %+5.3e %c %+2.0e%c  %2.0e%c  %2.0e%c  %2.0e%c  %2.0e%c  %6.4f   %d %d %d\n",(int)info->iter, 32,info->pcost, 32,info->dcost, 32, info->gap, 32, info->pres, 32, info->dres, 32, info->kapovert, 32, info->mu, 32, info->step, (int)info->nitref1, (int)info->nitref2, (int)info->nitref3);
#endif
	}

/* enable to flush printf in Matlab immediately */
#ifdef MATLAB_MEX_FILE
#if defined MATLAB_FLUSH_PRINTS
    mexEvalString("pause(0.0001);");
#endif
#endif
 
    
}

void deleteLastProgressLine( stats* info )
{
    idxint i;
    idxint offset = 0;
    
    if( info->kapovert < 0 ) offset++;
    if( info->mu < 0) offset++;
    if( info->pres < 0 ) offset++;
    if (info->dres < 0 ) offset++;
    
    for (i=0; i<82+offset; i++) {
        PRINTTEXT("%c",8);
    }
}
#endif



/**
 * Prepares the affine RHS for KKT system.
 * Given the special way we store the KKT matrix (sparse representation
 * of the scalings for the second-order cone), we need this to prepare
 * the RHS before solving the KKT system in the special format.
 */
void RHS_affine(pwork* w)
{
	pfloat* RHS = w->KKT->RHS2; 
	idxint n = w->n;
	idxint p = w->p;
	idxint i, j, k, l;
	idxint* Pinv = w->KKT->Pinv;

	j = 0;
	for( i=0; i < n; i++ ){ RHS[Pinv[j++]] = w->rx[i]; }
	for( i=0; i < p; i++ ){ RHS[Pinv[j++]] = -w->ry[i]; }
	for( i=0; i < w->C->lpc->p; i++ ){ RHS[Pinv[j++]] = w->s[i] - w->rz[i]; }
	k = w->C->lpc->p;
	for( l=0; l < w->C->nsoc; l++ ){
		for( i=0; i < w->C->soc[l].p; i++ ){ 
			RHS[Pinv[j++]] = w->s[k] - w->rz[k]; k++;			
		}
#if CONEMODE == 0
		RHS[Pinv[j++]] = 0;
        RHS[Pinv[j++]] = 0;
#endif
	}	
}



/**
 * Prepares the RHS for computing the combined search direction.
 */
void RHS_combined(pwork* w)
{
	pfloat* ds1 = w->KKT->work1;
	pfloat* ds2 = w->KKT->work2;
	idxint i, j, k, l;
	pfloat sigmamu = w->info->sigma * w->info->mu;
	pfloat one_minus_sigma = 1.0 - w->info->sigma;
	idxint* Pinv = w->KKT->Pinv;

	/* ds = lambda o lambda + W\s o Wz - sigma*mu*e) */	
	conicProduct(w->lambda, w->lambda, w->C, ds1);
	conicProduct(w->dsaff_by_W, w->W_times_dzaff, w->C, ds2);	
	for( i=0; i < w->C->lpc->p; i++ ){ ds1[i] += ds2[i] - sigmamu; }
	k = w->C->lpc->p;
	for( i=0; i < w->C->nsoc; i++ ){
		ds1[k] += ds2[k] - sigmamu; k++;
		for( j=1; j < w->C->soc[i].p; j++ ){ ds1[k] += ds2[k]; k++; }
	}

	/* dz = -(1-sigma)*rz + W*(lambda \ ds) */
	conicDivision(w->lambda, ds1, w->C, w->dsaff_by_W);
	scale(w->dsaff_by_W, w->C, ds1);
	
	/* copy in RHS */
	j = 0;
	for( i=0; i < w->n; i++ ){ w->KKT->RHS2[Pinv[j++]] *= one_minus_sigma; }
	for( i=0; i < w->p; i++ ){ w->KKT->RHS2[Pinv[j++]] *= one_minus_sigma; }
    for( i=0; i < w->C->lpc->p; i++) { w->KKT->RHS2[Pinv[j++]] = -one_minus_sigma*w->rz[i] + ds1[i]; }
	k = w->C->lpc->p;
	for( l=0; l < w->C->nsoc; l++ ){
		for( i=0; i < w->C->soc[l].p; i++ ){ 
			w->KKT->RHS2[Pinv[j++]] = -one_minus_sigma*w->rz[k] + ds1[k];
			k++;
		}
#if CONEMODE == 0
		w->KKT->RHS2[Pinv[j++]] = 0;
        w->KKT->RHS2[Pinv[j++]] = 0;
#endif
	}	
}



/*
 * Line search according to Vandenberghe (cf. �8 in his manual).
 */
pfloat lineSearch(pfloat* lambda, pfloat* ds, pfloat* dz, pfloat tau, pfloat dtau, pfloat kap, pfloat dkap, cone* C, kkt* KKT)
{
	idxint i, j, cone_start, conesize;
	pfloat rhomin, sigmamin, alpha, lknorm, lknorminv, rhonorm, sigmanorm, conic_step, temp;
	pfloat lkbar_times_dsk, lkbar_times_dzk, factor;
	pfloat* lk;
	pfloat* dsk;
	pfloat* dzk;
	pfloat* lkbar = KKT->work1;
	pfloat* rho = KKT->work2;
	pfloat* sigma = KKT->work2;
    pfloat minus_tau_by_dtau = -tau/dtau;
    pfloat minus_kap_by_dkap = -kap/dkap;
    

	/* LP cone */    
	if( C->lpc->p > 0 ){
		rhomin = ds[0] / lambda[0];  sigmamin = dz[0] / lambda[0];
		for( i=1; i < C->lpc->p; i++ ){
			rho[0] = ds[i] / lambda[i];   if( rho[0] < rhomin ){ rhomin = rho[0]; }
			sigma[0] = dz[i] / lambda[i]; if( sigma[0] < sigmamin ){ sigmamin = sigma[0]; }
		}

		if( -sigmamin > -rhomin ){ 
			alpha = sigmamin < 0 ? 1.0 / (-sigmamin) : 1.0 / EPS;
		} else {
			alpha = rhomin < 0 ? 1.0 / (-rhomin) : 1.0 / EPS;
		}
    } else {
        alpha = 10;
    }

    /* tau and kappa */
    if( minus_tau_by_dtau > 0 && minus_tau_by_dtau < alpha )
    {
        alpha = minus_tau_by_dtau;
    }
    if( minus_kap_by_dkap > 0 && minus_kap_by_dkap < alpha )
    {
        alpha = minus_kap_by_dkap;
    }
        

	/* Second-order cone */
	cone_start = C->lpc->p;
	for( i=0; i < C->nsoc; i++ ){
		
		/* indices */
		conesize = C->soc[i].p;
		lk = lambda + cone_start;  dsk = ds + cone_start;  dzk = dz + cone_start;
		
		/* normalize */
		lknorm = sqrt( lk[0]*lk[0] - eddot(conesize-1, lk+1, lk+1) );
		for( j=0; j < conesize; j++ ){ lkbar[j] = lk[j] / lknorm; }
		lknorminv = 1.0 / lknorm;

		/* calculate products */
		lkbar_times_dsk = lkbar[0]*dsk[0];
		for( j=1; j < conesize; j++ ){ lkbar_times_dsk -= lkbar[j]*dsk[j]; }
		lkbar_times_dzk = lkbar[0]*dzk[0];
		for( j=1; j < conesize; j++ ){ lkbar_times_dzk -= lkbar[j]*dzk[j]; }

		/* now construct rhok and sigmak, the first element is different */
		rho[0] = lknorminv * lkbar_times_dsk;
		factor = (lkbar_times_dsk+dsk[0])/(lkbar[0]+1);
		for( j=1; j < conesize; j++ ){ rho[j] = lknorminv*(dsk[j] - factor*lkbar[j]); }
		rhonorm = norm2(rho+1, conesize-1) - rho[0];

		sigma[0] = lknorminv * lkbar_times_dzk;
		factor = (lkbar_times_dzk+dzk[0])/(lkbar[0]+1);
		for( j=1; j < conesize; j++ ){ sigma[j] = lknorminv*(dzk[j] - factor*lkbar[j]); }
		sigmanorm = norm2(sigma+1, conesize-1) - sigma[0];
		
		/* update alpha */
		conic_step = 0;
		if( rhonorm > conic_step ){ conic_step = rhonorm; }
		if( sigmanorm > conic_step ){ conic_step = sigmanorm; }
		if( conic_step != 0 ){
			temp = 1.0 / conic_step;
			if( temp < alpha ){ alpha = temp; }
		}

		cone_start += C->soc[i].p;
		
	}	
    
    /* saturate between STEPMIN and STEPMAX */
    if( alpha > STEPMAX ) alpha = STEPMAX;
    if( alpha < STEPMIN ) alpha = STEPMIN;
    
    /* return alpha */
	return alpha;
}



/**
 * Scales variables by 1.0/tau, i.e. computes 
 * x = x./tau, y = y./tau, z = z./tau, s = s./tau
 */
void backscale(pwork *w)
{
	idxint i;
#if defined EQUILIBRATE && EQUILIBRATE > 0
    /* We performed a change of variables on x so this unsets the equilibration */
	for( i=0; i < w->n; i++ ){ w->x[i] /= (w->xequil[i] * w->tau); }
    for( i=0; i < w->p; i++ ){ w->y[i] /= (w->Aequil[i] * w->tau); }
	for( i=0; i < w->m; i++ ){ w->z[i] /= (w->Gequil[i] * w->tau); }
	for( i=0; i < w->m; i++ ){ w->s[i] /= (w->Gequil[i] * w->tau); }
#else
    /* standard back scaling without equilibration */
    for( i=0; i < w->n; i++ ){ w->x[i] /= w->tau; }
    for( i=0; i < w->p; i++ ){ w->y[i] /= w->tau; }
	for( i=0; i < w->m; i++ ){ w->z[i] /= w->tau; }
	for( i=0; i < w->m; i++ ){ w->s[i] /= w->tau; }
#endif
}



/*
 * Main solver routine.
 */
idxint ECOS_solve(pwork* w)
{
	idxint i, initcode, KKT_FACTOR_RETURN_CODE;
	pfloat dtau_denom, dtauaff, dkapaff, sigma, dtau, dkap, bkap, pres_prev;
	idxint exitcode = ECOS_FATAL;
    
#if DEBUG
    char fn[20];
#endif

#if (defined _WIN32 || defined _WIN64 )
	/* sets width of exponent for floating point numbers to 2 instead of 3 */
	unsigned int old_output_format = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

#if PROFILING > 0
	timer tsolve;
#endif
#if PROFILING > 1
    timer tfactor, tkktsolve;
#endif
    
#if PROFILING > 0
    /* start timer */
    tic(&tsolve);
#endif
	
	/* Initialize solver */
    initcode = init(w);
	if( initcode == ECOS_FATAL ){
#if PRINTLEVEL > 0
        if( w->stgs->verbose ) PRINTTEXT("\nFatal error during initialization, aborting.");
#endif
        return ECOS_FATAL;
    }
    
    
    
	/* MAIN INTERIOR POINT LOOP ---------------------------------------------------------------------- */
	for( w->info->iter = 0; w->info->iter <= w->stgs->maxit; w->info->iter++ ){
        
		/* Compute residuals */
		computeResiduals(w);
        
		/* Update statistics */
		updateStatistics(w);

#if PRINTLEVEL > 1
		/* Print info */
		if( w->stgs->verbose ) printProgress(w->info);
#endif
        
        /* SAFEGUARD: Backtrack to best previously seen iterate if
         *
         * - the update was bad such that the primal residual PRES has increased by a factor of SAFEGUARD, or
         * - the gap became negative
         *
         * If the safeguard is activated, the solver tests if reduced precision has been reached, and reports
         * accordingly. If not even reduced precision is reached, ECOS returns the flag ECOS_NUMERICS.
         */
        if( w->info->iter > 0 && (w->info->pres > SAFEGUARD*pres_prev || w->info->gap < 0) ){
#if PRINTLEVEL > 1
            if( w->stgs->verbose ) deleteLastProgressLine( w->info );
            if( w->stgs->verbose ) PRINTTEXT("Unreliable search direction detected, recovering best iterate (%d) and stopping.\n", (int)w->best_info->iter);
#endif
            restoreBestIterate( w );
            
            /* Determine whether we have reached at least reduced accuracy */
            exitcode = checkExitConditions( w, ECOS_INACC_OFFSET );
            
            /* if not, exit anyways */
            if( exitcode == ECOS_NOT_CONVERGED_YET ){
                exitcode = ECOS_NUMERICS;
#if PRINTLEVEL > 0
                if( w->stgs->verbose ) PRINTTEXT("\nNUMERICAL PROBLEMS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", MAX(w->info->dres, w->info->pres), w->info->relgap, w->info->gap);
#endif
                break;
            } else {
                break;
            }
        }
        pres_prev = w->info->pres;
        

		/* Check termination criteria to full precision and exit if necessary */
		exitcode = checkExitConditions( w, 0 );
        if( exitcode == ECOS_NOT_CONVERGED_YET ){
            
            /*
             * Full precision has not been reached yet. Check for two more cases of exit:
             *  (i) min step size, in which case we assume we won't make progress any more, and
             * (ii) maximum number of iterations reached
             * If these two are not fulfilled, another iteration will be made.
             */
            
            /* Did the line search cock up? (zero step length) */
            if( w->info->iter > 0 && w->info->step == STEPMIN*GAMMA ){
#if PRINTLEVEL > 0
                if( w->stgs->verbose ) deleteLastProgressLine( w->info );
                if( w->stgs->verbose ) PRINTTEXT("No further progress possible, recovering best iterate (%d) and stopping.", (int)w->best_info->iter );
#endif
                restoreBestIterate( w );
                
                /* Determine whether we have reached reduced precision */
                exitcode = checkExitConditions( w, ECOS_INACC_OFFSET );
                if( exitcode == ECOS_NOT_CONVERGED_YET ){
                    exitcode = ECOS_NUMERICS;
#if PRINTLEVEL > 0
                    if( w->stgs->verbose ) PRINTTEXT("\nNUMERICAL PROBLEMS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", MAX(w->info->dres, w->info->pres), w->info->relgap, w->info->gap);
#endif
                }
                break;
            }
            /* MAXIT reached? */
            else if( w->info->iter == w->stgs->maxit ){
                
                /* Determine whether current iterate is better than what we had so far */
                if( compareStatistics( w->info, w->best_info) ){
#if PRINTLEVEL > 0
                    if( w->stgs->verbose ) PRINTTEXT("Maximum number of iterations reached, stopping.\n");
#endif
                } else
                {
#if PRINTLEVEL > 0
                    if( w->stgs->verbose ) PRINTTEXT("Maximum number of iterations reached, recovering best iterate (%d) and stopping.\n", (int)w->best_info->iter);
#endif
                    restoreBestIterate( w );
                }
                
                /* Determine whether we have reached reduced precision */
                exitcode = checkExitConditions( w, ECOS_INACC_OFFSET );
                if( exitcode == ECOS_NOT_CONVERGED_YET ){
                    exitcode = ECOS_MAXIT;
#if PRINTLEVEL > 0
                    if( w->stgs->verbose ) PRINTTEXT("\nRAN OUT OF ITERATIONS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", MAX(w->info->dres, w->info->pres), w->info->relgap, w->info->gap);
#endif
                }
                break;
            }
        } else {
            
            /* Full precision has been reached, stop solver */
            break;
        }
        
		
        
        /* SAFEGUARD:
         * Check whether current iterate is worth keeping as the best solution so far,
         * before doing another iteration
         */
        if (w->info->iter == 0) {
            /* we're at the first iterate, so there's nothing to compare yet */
            saveIterateAsBest( w );
        } else if( compareStatistics( w->info, w->best_info) ){
            /* PRINTTEXT("Better solution found, saving as best so far \n"); */
            saveIterateAsBest( w );
        }
        

		/* Compute scalings */
		if( updateScalings(w->C, w->s, w->z, w->lambda) == OUTSIDE_CONE ){
            
            /* SAFEGUARD: we have to recover here */
#if PRINTLEVEL > 0
            if( w->stgs->verbose ) deleteLastProgressLine( w->info );
            if( w->stgs->verbose ) PRINTTEXT("Slacks/multipliers leaving the cone, recovering best iterate (%d) and stopping.\n", (int)w->best_info->iter);
#endif
            restoreBestIterate( w );
            
            /* Determine whether we have reached at least reduced accuracy */
            exitcode = checkExitConditions( w, ECOS_INACC_OFFSET );
            if( exitcode == ECOS_NOT_CONVERGED_YET ){
#if PRINTLEVEL > 0
                if( w->stgs->verbose ) PRINTTEXT("\nNUMERICAL PROBLEMS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", MAX(w->info->dres, w->info->pres), w->info->relgap, w->info->gap);
#endif
                return ECOS_OUTCONE;

            } else {
                break;
            }
        }
        
		/* Update KKT matrix with scalings */
		kkt_update(w->KKT->PKPt, w->KKT->PK, w->C);
        
#if DEBUG > 0
        /* DEBUG: Store matrix to be factored */
        sprintf(fn, "PKPt_updated_%02i.txt", (int)w->info->iter);
        dumpSparseMatrix(w->KKT->PKPt, fn);
#endif
        /* factor KKT matrix */
#if PROFILING > 1
		tic(&tfactor);
        KKT_FACTOR_RETURN_CODE = kkt_factor(w->KKT, w->stgs->eps, w->stgs->delta, &w->info->tfactor_t1, &w->info->tfactor_t2);
        w->info->tfactor += toc(&tfactor);
#else
        KKT_FACTOR_RETURN_CODE = kkt_factor(w->KKT, w->stgs->eps, w->stgs->delta);
#endif
        
#if DEBUG > 0
        /* DEBUG: store factor */
        sprintf(fn, "PKPt_factor_%02i.txt", (int)w->info->iter);
        dumpSparseMatrix(w->KKT->L, fn);
#endif

		/* Solve for RHS1, which is used later also in combined direction */
#if PROFILING > 1
		tic(&tkktsolve);
#endif
		w->info->nitref1 = kkt_solve(w->KKT, w->A, w->G, w->KKT->RHS1, w->KKT->dx1, w->KKT->dy1, w->KKT->dz1, w->n, w->p, w->m, w->C, 0, w->stgs->nitref);
#if PROFILING > 1
		w->info->tkktsolve += toc(&tkktsolve);
#endif
        
#if DEBUG > 0 && PRINTLEVEL > 2
        /* Print result of linear system solve */
        printDenseMatrix(w->KKT->dx1, 1, 5, "dx1(1:5)");
        printDenseMatrix(w->KKT->dy1, 1, 5, "dy1(1:5)");
        printDenseMatrix(w->KKT->dz1, 1, 5, "dz1(1:5)");
#endif
  
		/* AFFINE SEARCH DIRECTION (predictor, need dsaff and dzaff only) */
		RHS_affine(w);
#if PROFILING > 1
		tic(&tkktsolve);
#endif
		w->info->nitref2 = kkt_solve(w->KKT, w->A, w->G, w->KKT->RHS2, w->KKT->dx2, w->KKT->dy2, w->KKT->dz2, w->n, w->p, w->m, w->C, 0, w->stgs->nitref);
#if PROFILING > 1
		w->info->tkktsolve += toc(&tkktsolve);
#endif
        
		/* dtau_denom = kap/tau - (c'*x1 + by1 + h'*z1); */
		dtau_denom = w->kap/w->tau - eddot(w->n, w->c, w->KKT->dx1) - eddot(w->p, w->b, w->KKT->dy1) - eddot(w->m, w->h, w->KKT->dz1);
		
        /* dtauaff = (dt + c'*x2 + by2 + h'*z2) / dtau_denom; */
		dtauaff = (w->rt - w->kap + eddot(w->n, w->c, w->KKT->dx2) + eddot(w->p, w->b, w->KKT->dy2) + eddot(w->m, w->h, w->KKT->dz2)) / dtau_denom;
        
		/* dzaff = dz2 + dtau_aff*dz1 */
		for( i=0; i<w->m; i++ ){ w->W_times_dzaff[i] = w->KKT->dz2[i] + dtauaff*w->KKT->dz1[i]; } 
		scale(w->W_times_dzaff, w->C, w->W_times_dzaff);

		/* W\dsaff = -W*dzaff -lambda; */		
		for( i=0; i<w->m; i++ ){ w->dsaff_by_W[i] = -w->W_times_dzaff[i] - w->lambda[i]; }
		
		/* dkapaff = -(bkap + kap*dtauaff)/tau; bkap = kap*tau*/
		dkapaff = -w->kap - w->kap/w->tau*dtauaff;
        
        /* Line search on W\dsaff and W*dzaff */
		w->info->step_aff = lineSearch(w->lambda, w->dsaff_by_W, w->W_times_dzaff, w->tau, dtauaff, w->kap, dkapaff, w->C, w->KKT);
        
		/* Centering parameter */
        sigma = 1.0 - w->info->step_aff;
        sigma = sigma*sigma*sigma;
        if( sigma > SIGMAMAX ) sigma = SIGMAMAX;
        if( sigma < SIGMAMIN ) sigma = SIGMAMIN;
        w->info->sigma = sigma;
        
		
		/* COMBINED SEARCH DIRECTION */
		RHS_combined(w);
#if PROFILING > 1
		tic(&tkktsolve);
#endif
		w->info->nitref3 = kkt_solve(w->KKT, w->A, w->G, w->KKT->RHS2, w->KKT->dx2, w->KKT->dy2, w->KKT->dz2, w->n, w->p, w->m, w->C, 0, w->stgs->nitref);
#if PROFILING > 1
		w->info->tkktsolve += toc(&tkktsolve);
#endif
        
  		/* bkap = kap*tau + dkapaff*dtauaff - sigma*info.mu; */
		bkap = w->kap*w->tau + dkapaff*dtauaff - sigma*w->info->mu;

		/* dtau = ((1-sigma)*rt - bkap/tau + c'*x2 + by2 + h'*z2) / dtau_denom; */		
		dtau = ((1-sigma)*w->rt - bkap/w->tau + eddot(w->n, w->c, w->KKT->dx2) + eddot(w->p, w->b, w->KKT->dy2) + eddot(w->m, w->h, w->KKT->dz2)) / dtau_denom;
      	
		/* dx = x2 + dtau*x1;     dy = y2 + dtau*y1;       dz = z2 + dtau*z1; */
		for( i=0; i < w->n; i++ ){ w->KKT->dx2[i] += dtau*w->KKT->dx1[i]; }
		for( i=0; i < w->p; i++ ){ w->KKT->dy2[i] += dtau*w->KKT->dy1[i]; }
		for( i=0; i < w->m; i++ ){ w->KKT->dz2[i] += dtau*w->KKT->dz1[i]; }

		/*  ds_by_W = -(lambda \ bs + conelp_timesW(scaling,dz,dims)); */
		/* note that ath this point w->dsaff_by_W holds already (lambda \ ds) */
		scale(w->KKT->dz2, w->C, w->W_times_dzaff);
		for( i=0; i < w->m; i++ ){ w->dsaff_by_W[i] = -(w->dsaff_by_W[i] + w->W_times_dzaff[i]); }

		/* dkap = -(bkap + kap*dtau)/tau; */
		dkap = -(bkap + w->kap*dtau)/w->tau;

		/* Line search on combined direction */
		w->info->step = lineSearch(w->lambda, w->dsaff_by_W, w->W_times_dzaff, w->tau, dtau, w->kap, dkap, w->C, w->KKT) * w->stgs->gamma;
		
		/* ds = W*ds_by_W */
		scale(w->dsaff_by_W, w->C, w->dsaff);

		/* Update variables */
		for( i=0; i < w->n; i++ ){ w->x[i] += w->info->step * w->KKT->dx2[i]; }
		for( i=0; i < w->p; i++ ){ w->y[i] += w->info->step * w->KKT->dy2[i]; }
		for( i=0; i < w->m; i++ ){ w->z[i] += w->info->step * w->KKT->dz2[i]; }
		for( i=0; i < w->m; i++ ){ w->s[i] += w->info->step * w->dsaff[i]; }
		w->kap += w->info->step * dkap;
		w->tau += w->info->step * dtau;
	}

	/* scale variables back */    
	backscale(w);

	/* stop timer */
#if PROFILING > 0
	w->info->tsolve = toc(&tsolve);
#endif

#if PRINTLEVEL > 0
#if PROFILING > 0
	if( w->stgs->verbose ) PRINTTEXT("\nRuntime: %f seconds.", w->info->tsetup + w->info->tsolve);
#endif
	if( w->stgs->verbose ) PRINTTEXT("\n\n");
#endif

	return exitcode;
}


/*
 * Updates one element of the RHS vector h of inequalities
 * After the call, w->h[idx] = value (but equilibrated)
 */
void ecos_updateDataEntry_h(pwork* w, idxint idx, pfloat value)
{
#if defined EQUILIBRATE && EQUILIBRATE > 0
    w->h[idx] = value / w->Gequil[idx];
#else
    w->h[idx] = value;
#endif
}
