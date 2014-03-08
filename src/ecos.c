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


/*
 * Initializes the solver.
 *
 * This function assumes that the KKT matrix K is already in the form
 *
 *		[0  A'  G']
 * K =  [A  0   0 ]
 *      [G  0  -I ]
 *
 * The preprocessor/codegen takes care of this, we just have to be aware
 * of this implicit assumption.
 */
idxint init(pwork* w)
{
	idxint i, KKT_FACTOR_RETURN_CODE;
	idxint* Pinv = w->KKT->Pinv;

#if PROFILING > 1
	timer tfactor, tkktsolve;
#endif

	w->KKT->delta = w->stgs->delta;

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
	w->cx = ddot(w->n, w->c, w->x);
	w->by = w->p > 0 ? ddot(w->p, w->b, w->y) : 0.0;
	w->hz = ddot(w->m, w->h, w->z);
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
	info->gap = ddot(w->m, w->s, w->z);
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
		PRINTTEXT("\nECOS %s - (c) A. Domahidi, Automatic Control Lab, ETH Zurich, 2012-2014.\n\n", ECOS_VERSION);
#endif
#if defined _WIN32 || defined _WIN64		
		PRINTTEXT("It     pcost       dcost      gap   pres   dres    k/t    mu     step    IR\n");
		PRINTTEXT("%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e   N/A    %d %d -\n",(int)info->iter, info->pcost, info->dcost, info->gap, info->pres, info->dres, info->kapovert, info->mu, (int)info->nitref1, (int)info->nitref2);
#else
		PRINTTEXT("It     pcost         dcost      gap     pres    dres     k/t     mu     step    IR\n");
		PRINTTEXT("%2d  %c%+5.3e  %c%+5.3e  %c%+2.0e  %c%2.0e  %c%2.0e  %c%2.0e  %c%2.0e   N/A    %d %d -\n",(int)info->iter, 32, info->pcost, 32, info->dcost, 32, info->gap, 32, info->pres, 32, info->dres, 32, info->kapovert, 32, info->mu, (int)info->nitref1, (int)info->nitref2);
#endif	
	}  else {
#if defined _WIN32 || defined _WIN64
		PRINTTEXT("%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e  %6.4f  %d %d %d\n",(int)info->iter, info->pcost, info->dcost, info->gap, info->pres, info->dres, info->kapovert, info->mu, info->step, (int)info->nitref1, (int)info->nitref2, (int)info->nitref3);
#else
		PRINTTEXT("%2d  %c%+5.3e%c  %+5.3e %c %+2.0e%c  %2.0e%c  %2.0e%c  %2.0e%c  %2.0e  %6.4f  %d %d %d\n",(int)info->iter, 32,info->pcost, 32,info->dcost, 32, info->gap, 32, info->pres, 32, info->dres, 32, info->kapovert, 32, info->mu, info->step, (int)info->nitref1, (int)info->nitref2, (int)info->nitref3);
#endif
	}

/* enable to flush printf in Matlab immediately */
#if PRINTLEVEL > 0
#ifdef MATLAB_MEX_FILE
#if defined MATLAB_FLUSH_PRINTS
    mexEvalString("drawnow;");
#endif
#endif
#endif
 
    
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
	for( i=0; i < p; i++ ){ RHS[Pinv[j++]] = w->ry[i]; }
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
 * Line search according to Vandenberghe (cf. §8 in his manual).
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
		lknorm = sqrt( lk[0]*lk[0] - ddot(conesize-1, lk+1, lk+1) );
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
    /* We performed a change of variables on x so this recovers it also */
	for( i=0; i < w->n; i++ ){ w->x[i] /= (w->xequil[i] * w->tau); }
	for( i=0; i < w->p; i++ ){ w->y[i] /= (w->Aequil[i] * w->tau); }
	for( i=0; i < w->m; i++ ){ w->z[i] /= (w->Gequil[i] * w->tau); }
	for( i=0; i < w->m; i++ ){ w->s[i] /= (w->Gequil[i] * w->tau); }
}



/*
 * Main solver routine.
 */
idxint ECOS_solve(pwork* w)
{
	idxint i, initcode, KKT_FACTOR_RETURN_CODE;
	pfloat dtau_denom, dtauaff, dkapaff, sigma, dtau, dkap, bkap, pres_prev;
	idxint exitcode = ECOS_FATAL;

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
        
        
        /* SAFEGUARD: Backtrack to old iterate if the update was bad such that the primal residual PRES has
         *            increased by a factor of SAFEGUARD.
         * If the safeguard is activated, the solver quits with the flag ECOS_NUMERICS.
         */
        if( w->info->iter > 0 && w->info->pres > SAFEGUARD*pres_prev ){
#if PRINTLEVEL > 1
            if( w->stgs->verbose ) PRINTTEXT("\nNUMERICAL PROBLEMS, recovering iterate %d and stopping.\n", (int)w->info->iter-1);
#endif
            /* Backtrack */
            for( i=0; i < w->n; i++ ){ w->x[i] -= w->info->step * w->KKT->dx2[i]; }
            for( i=0; i < w->p; i++ ){ w->y[i] -= w->info->step * w->KKT->dy2[i]; }
            for( i=0; i < w->m; i++ ){ w->z[i] -= w->info->step * w->KKT->dz2[i]; }
            for( i=0; i < w->m; i++ ){ w->s[i] -= w->info->step * w->dsaff[i]; }
            w->kap -= w->info->step * dkap;
            w->tau -= w->info->step * dtau;
            exitcode = ECOS_NUMERICS;
            computeResiduals(w);
            updateStatistics(w);
            break;
        }
        pres_prev = w->info->pres;
        

		/* Check termination criteria and exit if necessary */
		/* Optimal? */
		if( ( ( -w->cx > 0 ) || ( -w->by - w->hz > 0) ) &&
            w->info->pres < w->stgs->feastol && w->info->dres < w->stgs->feastol &&
			( w->info->gap < w->stgs->abstol || w->info->relgap < w->stgs->reltol ) ){
#if PRINTLEVEL > 0
			if( w->stgs->verbose ) PRINTTEXT("\nOPTIMAL (within feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", w->stgs->feastol, w->stgs->reltol, w->stgs->abstol);
#endif
	        exitcode = ECOS_OPTIMAL;
			break;
		}            
		/* Primal infeasible? */
		else if( ((w->info->pinfres != NAN) && (w->info->pinfres < w->stgs->feastol)) ||
                 ((w->tau < w->stgs->feastol) && (w->kap < w->stgs->feastol && w->info->pinfres < w->stgs->feastol)) ){
#if PRINTLEVEL > 0
			if( w->stgs->verbose ) PRINTTEXT("\nPRIMAL INFEASIBLE (within feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", w->stgs->feastol, w->stgs->reltol, w->stgs->abstol);
#endif
			w->info->pinf = 1;
			w->info->dinf = 0;
			exitcode = ECOS_PINF;
			break;
		}        
		/* Dual infeasible? */
		else if( (w->info->dinfres != NAN) && (w->info->dinfres < w->stgs->feastol) ){
#if PRINTLEVEL > 0
			if( w->stgs->verbose ) PRINTTEXT("\nUNBOUNDED (within feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).", w->stgs->feastol, w->stgs->reltol, w->stgs->abstol);
#endif
			w->info->pinf = 0;  
			w->info->dinf = 1;        
			exitcode = ECOS_DINF;
			break;
		}   
		/* Did the line search cock up? (zero step length) */
		else if( w->info->iter > 0 && w->info->step == STEPMIN*GAMMA ){
#if PRINTLEVEL > 0
			if( w->stgs->verbose ) PRINTTEXT("\nNo further progress possible (- numerics?), exiting.");
#endif
			exitcode = ECOS_NUMERICS;
			break;
		}
		/* MAXIT reached? */
		else if( w->info->iter == w->stgs->maxit ){
#if PRINTLEVEL > 0
			if( w->stgs->verbose ) PRINTTEXT("\nMaximum number of iterations reached, exiting.");
#endif
			exitcode = ECOS_MAXIT;
			break;
		}


		/* Compute scalings */
		if( updateScalings(w->C, w->s, w->z, w->lambda) == OUTSIDE_CONE ){
#if PRINTLEVEL > 0
            if( w->stgs->verbose ) PRINTTEXT("\nSlacks or multipliers leaving the positive orthant (- numerics ?), exiting.\n");
#endif
            return ECOS_OUTCONE;
        }
        
		/* Update KKT matrix with scalings */
		kkt_update(w->KKT->PKPt, w->KKT->PK, w->C);
        
#if DEBUG > 0
        dumpSparseMatrix(w->KKT->PKPt,"PKPt_updated.txt");
#endif
        /* factor KKT matrix */
#if PROFILING > 1
		tic(&tfactor);
        KKT_FACTOR_RETURN_CODE = kkt_factor(w->KKT, w->stgs->eps, w->stgs->delta, &w->info->tfactor_t1, &w->info->tfactor_t2);
        w->info->tfactor += toc(&tfactor);
#else
        KKT_FACTOR_RETURN_CODE = kkt_factor(w->KKT, w->stgs->eps, w->stgs->delta);
#endif

		/* Solve for RHS1, which is used later also in combined direction */
#if PROFILING > 1
		tic(&tkktsolve);
#endif
		w->info->nitref1 = kkt_solve(w->KKT, w->A, w->G, w->KKT->RHS1, w->KKT->dx1, w->KKT->dy1, w->KKT->dz1, w->n, w->p, w->m, w->C, 0, w->stgs->nitref);
#if PROFILING > 1
		w->info->tkktsolve += toc(&tkktsolve);
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
		dtau_denom = w->kap/w->tau - ddot(w->n, w->c, w->KKT->dx1) - ddot(w->p, w->b, w->KKT->dy1) - ddot(w->m, w->h, w->KKT->dz1);
		
        /* dtauaff = (dt + c'*x2 + by2 + h'*z2) / dtau_denom; */
		dtauaff = (w->rt - w->kap + ddot(w->n, w->c, w->KKT->dx2) + ddot(w->p, w->b, w->KKT->dy2) + ddot(w->m, w->h, w->KKT->dz2)) / dtau_denom;
        
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
		dtau = ((1-sigma)*w->rt - bkap/w->tau + ddot(w->n, w->c, w->KKT->dx2) + ddot(w->p, w->b, w->KKT->dy2) + ddot(w->m, w->h, w->KKT->dz2)) / dtau_denom;
      	
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









