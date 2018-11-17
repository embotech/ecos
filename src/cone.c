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


/* cone functions */
#include "cone.h"
#include "spla.h"
#include "ecos.h"
#include "expcone.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* PRIVATE METHODS ===================================================== */
/*
 * Computes u0^2 - u1'*u1 for second order cones
 */
pfloat socres(pfloat* u, idxint p)
{
	pfloat res = u[0]*u[0];
	idxint i;
	for( i=1; i<p; i++) { res -= u[i]*u[i]; }
	return res;
}


/* PUBLIC METHODS ====================================================== */

/**
 * Scales a conic variable such that it lies strictly in the cone.
 * If it is already in the cone, r is simply copied to s.
 * Otherwise s = r + (1+alpha)*e where alpha is the biggest residual.
 */
void bring2cone(cone* C, pfloat* r, pfloat* s)
{
	pfloat alpha = -GAMMA;
	pfloat cres, r1square;
	idxint i, l, j;

	/*
	 * 1. Find maximum residual ----------------------------------------------
	 */

	/* LP cone */
	for( i=0; i<C->lpc->p; i++ ){
		if( r[i] <= 0 && -r[i] > alpha ){ alpha = -r[i]; }
	}

	/* Second-order cone */
	for( l=0; l < C->nsoc; l++ ){
		cres = r[i++]; r1square = 0;
		for( j=1; j<C->soc[l].p; j++ ){ r1square += r[i]*r[i]; i++; }
		cres -= sqrt(r1square);
		if( cres <= 0 && -cres > alpha ){ alpha = -cres; }
	}


	/*
	 * 2. compute s = r + (1+alpha)*e -------------------------------------
	 */

	alpha += 1.0;

	/* LP cone */
	for( i=0; i < C->lpc->p; i++ ){
		s[i] = r[i] + alpha;
	}

	/* Second-order cone */
	for( l=0; l < C->nsoc; l++ ){
		s[i] = r[i] + alpha; i++;
		for( j=1; j < C->soc[l].p; j++ ){ s[i] = r[i]; i++; }
	}
}

#ifdef EXPCONE

/* Sets the initial point to the jordan algebra identity e times scaling for the symmetric cones
 * and the central ray for the exponential cone, scaled by scaling
 */
void unitInitialization(cone* C, pfloat* s, pfloat* z, pfloat scaling)
{
    idxint i,l,j;

    /* LP cone */
    for( i=0; i < C->lpc->p; i++ ){
        s[i] = scaling;
        z[i] = scaling;
    }

    /* Second-order cone */
    for( l=0; l < C->nsoc; l++ ){
        s[i] = scaling; z[i] = scaling;  i++;
        for( j=1; j < C->soc[l].p; j++ ){ s[i] = 0.0; z[i] = 0.0; i++; }
    }
    /* Exponential cone */
    for(l=0;l<C->nexc;l++)
    {
       s[i]   = scaling*(-1.051383945322714);
       s[i+1] = scaling*(1.258967884768947);
       s[i+2] = scaling*(0.556409619469370);
       z[i]   = scaling*(-1.051383945322714);
       z[i+1] = scaling*(1.258967884768947);
       z[i+2] = scaling*(0.556409619469370);
       i=i+3;
    }
}

#endif


/**
 * Update scalings.
 * Returns OUTSIDE_CONE as soon as any multiplier or slack leaves the cone,
 * as this indicates severe problems.
 */
#ifdef EXPCONE
idxint updateScalings(cone* C, pfloat* s, pfloat* z, pfloat* lambda, pfloat mu)
#else
idxint updateScalings(cone* C, pfloat* s, pfloat* z, pfloat* lambda)
#endif
{
	idxint i, l, k, p; /*, pm1; */
	pfloat sres, zres, snorm, znorm, gamma, one_over_2gamma;
	pfloat* sk;
	pfloat* zk;
    pfloat a, c, d, w, temp, divisor; /*, b; */
#if CONEMODE == 0
    pfloat u0, u0_square, u1, v1, d1, c2byu02_d, c2byu02, c_square;
#endif

	/* LP cone */
	for( i=0; i < C->lpc->p; i++ ){
		C->lpc->v[i] = SAFEDIV_POS(s[i], z[i]);
		C->lpc->w[i] = sqrt(C->lpc->v[i]);
	}

	/* Second-order cone */
	k = C->lpc->p;
	for( l=0; l < C->nsoc; l++ ){

		/* indices and variables */
		sk = s+k; zk = z+k; p = C->soc[l].p; /* pm1 = p-1; */

		/* check residuals and quit if they're negative */
		sres = socres(sk, p);  zres = socres(zk, p);
        if( sres <= 0 || zres <= 0 ){ return OUTSIDE_CONE; }

		/* normalize variables */
		snorm = sqrt(sres);    znorm = sqrt(zres);
		for( i=0; i<p; i++ ){ C->soc[l].skbar[i] = SAFEDIV_POS(sk[i],snorm); }
		for( i=0; i<p; i++ ){ C->soc[l].zkbar[i] = SAFEDIV_POS(zk[i],znorm); }
		C->soc[l].eta_square = SAFEDIV_POS(snorm,znorm);
		C->soc[l].eta = sqrt(C->soc[l].eta_square);

		/* Normalized Nesterov-Todd scaling point */
		gamma = 1.0;
		for( i=0; i<p; i++){ gamma += C->soc[l].skbar[i]*C->soc[l].zkbar[i]; }
		gamma = sqrt(0.5*gamma);
		one_over_2gamma = SAFEDIV_POS(0.5,gamma);
		a = one_over_2gamma*(C->soc[l].skbar[0] + C->soc[l].zkbar[0]);
		w = 0;
		for( i=1; i<p; i++ ){
			C->soc[l].q[i-1] = one_over_2gamma*(C->soc[l].skbar[i] - C->soc[l].zkbar[i]);
			w += C->soc[l].q[i-1]*C->soc[l].q[i-1];
		}
        C->soc[l].w = w;
        C->soc[l].a = a;

		/* pre-compute variables needed for KKT matrix (kkt_update uses those) */
        temp = 1.0 + a;
        /* b = SAFEDIV_POS(1.0,temp); */
        c = 1.0 + a + SAFEDIV_POS(w,temp);
        divisor = temp*temp;
        d = 1 + SAFEDIV_POS(2,temp) + SAFEDIV_POS(w,divisor);
#if CONEMODE > 0
        C->soc[l].c = c;
        C->soc[l].d = d;
#endif
#if CONEMODE == 0
        c_square = c*c;
        divisor = 1.0 + w*d;
        d1 = 0.5*(a*a + w*(1.0 - SAFEDIV_POS(c_square,divisor) ));  if( d1 < 0 ){ d1 = 0; }
        u0_square = a*a + w - d1;
        u0 = sqrt(u0_square);
        c2byu02 = SAFEDIV_POS((c*c),u0_square);
        c2byu02_d = c2byu02 - d;
        if (c2byu02_d <= 0) { return OUTSIDE_CONE; }
        v1 = sqrt(c2byu02_d);
        u1 = sqrt(c2byu02);
        C->soc[l].d1 = d1;
        C->soc[l].u0 = u0;
        C->soc[l].u1 = u1;
        C->soc[l].v1 = v1;
#endif

		/* increase offset for next cone */
		k += C->soc[l].p;
	}
#ifdef EXPCONE
       /*Exponential cones*/
        k = C->fexv;
        for(l=0;l<C->nexc;l++)
        {
            evalExpHessian(z+k, C->expc[l].v, mu);
            evalExpGradient(z+k, C->expc[l].g);
            k+=3;
        }
#endif
	/* lambda = W*z */
	scale(z, C, lambda);

	return INSIDE_CONE;
}

#ifdef EXPCONE
/* Evaulates log(s) + log(z) + log(t) + log(k) + logbarriersocp - (D+1) */
pfloat evalSymmetricBarrierValue(pfloat* siter, pfloat *ziter, pfloat tauIter, pfloat kapIter, cone* C, pfloat D)
{
   pfloat barrier = 0.0;
   idxint j,k,l;
   pfloat normAccumS = 0.0;
   pfloat normAccumZ = 0.0;
   idxint socDim;
   /* Positive orthant barrier */
   for(k=0;k<C->lpc->p;k++)
        barrier -= (siter[k]<=0||ziter[k]<=0)? ECOS_INFINITY : (log(siter[k])+log(ziter[k]));

    barrier-= (tauIter<=0||kapIter<=0) ? ECOS_INFINITY : (log(tauIter)+log(kapIter));
    /* Socp cones */
    for(l=0;l<C->nsoc;l++)
    {
        socDim = C->soc[l].p;
        normAccumS = 0.0;
        normAccumZ = 0.0;
        normAccumS = siter[k]*siter[k]; /* Root variable of the socp cone */
        normAccumZ = ziter[k]*ziter[k];
        k++;
        for(j=1;j<socDim;j++)
        {
            normAccumS -= siter[k]*siter[k];
            normAccumZ -= ziter[k]*ziter[k];
            k++;
        }
        barrier-= normAccumS<=0.0? ECOS_INFINITY : 0.5*log(normAccumS);
        barrier-= normAccumZ<=0.0? ECOS_INFINITY : 0.5*log(normAccumZ);

    }
    return barrier-D-1;
}
#endif

/**
 * Fast multiplication by scaling matrix.
 * Returns lambda = W*z
 * The exponential variables are not touched.
 */
void scale(pfloat* z, cone* C, pfloat* lambda)
{
	idxint i, j, l, cone_start;
	pfloat zeta, factor;

	/* LP cone */
	for( i=0; i < C->lpc->p; i++ ){ lambda[i] = C->lpc->w[i] * z[i]; }

	/* Second-order cone */
	cone_start = C->lpc->p;
	for( l=0; l < C->nsoc; l++ ){

		/* zeta = q'*z1 */
		zeta = 0;
		for( i=1; i < C->soc[l].p; i++ ){ zeta += C->soc[l].q[i-1] * z[cone_start + i]; }

		/* factor = z0 + zeta / (1+a); */
		factor = z[cone_start] + SAFEDIV_POS(zeta,(1+C->soc[l].a));

		/* second pass (on k): write out result */
		lambda[cone_start] = C->soc[l].eta*(C->soc[l].a*z[cone_start] + zeta); /* lambda[0] */
		for( i=1; i < C->soc[l].p; i++ ){
			j = cone_start+i;
			lambda[j] = C->soc[l].eta*(z[j] + factor*C->soc[l].q[i-1]);
		}

		cone_start += C->soc[l].p;
	}
}


/**
 *                                       [ D   v   u  ]
 * Fast multiplication with V =  eta^2 * [ v'  1   0  ] = W^2
 *                                       [ u   0  -1  ]
 * Computes y += W^2*x;
 */
void scale2add(pfloat *x, pfloat* y, cone* C)
{
    idxint i, l, cone_start, conesize, conesize_m1;
    pfloat *x1, *x2, *y1, *y2, eta_square, *q;
#if CONEMODE == 0
    pfloat *x3, *x4, *y3, *y4;
    pfloat d1, u0, u1, v1;
    pfloat v1x3_plus_u1x4;
    pfloat qtx2;
#else
    pfloat zeta, temp, a, w, c, d;
#endif

    /* LP cone */
	for( i=0; i < C->lpc->p; i++ ){ y[i] += C->lpc->v[i] * x[i]; }

    /* Second-order cone */
    cone_start = C->lpc->p;
#if CONEMODE == 0

	for( l=0; l < C->nsoc; l++ ){

        getSOCDetails(&C->soc[l], &conesize, &eta_square, &d1, &u0, &u1, &v1, &q);
        conesize_m1 = conesize - 1;

        x1 = x + cone_start;
        x2 = x1 + 1;
        x3 = x2 + conesize - 1;
        x4 = x3 + 1;

        y1 = y + cone_start;
        y2 = y1 + 1;
        y3 = y2 + conesize - 1;
        y4 = y3 + 1;

        /* y1 += d1*x1 + u0*x4 */
        y1[0] += eta_square*(d1*x1[0] + u0*x4[0]);

        /* y2 += x2 + v1*q*x3 + u1*q*x4 */
        v1x3_plus_u1x4 = v1*x3[0] + u1*x4[0];
        qtx2 = 0;
        for (i=0; i<conesize_m1; i++) {
            y2[i] += eta_square*(x2[i] + v1x3_plus_u1x4*q[i]);
            qtx2 += q[i]*x2[i];
        }

        /* y3 += v1*q'*x2 + x3 */
        y3[0] += eta_square*(v1*qtx2 + x3[0]);

        /* y4 += u0*x1 + u1*q'*x2 - x4 */
        y4[0] += eta_square*(u0*x1[0] + u1*qtx2 - x4[0]);

        /* prepare index for next cone */
        cone_start += conesize + 2;
    }

#else

	for( l=0; l < C->nsoc; l++ ){

        conesize = C->soc[l].p;
        conesize_m1 = conesize - 1;
        eta_square = C->soc[l].eta_square;
        a = C->soc[l].a;
        c = C->soc[l].c;
        d = C->soc[l].d;
        q = C->soc[l].q;
        w = C->soc[l].w;

        x1 = x + cone_start;
        x2 = x1 + 1;

        y1 = y + cone_start;
        y2 = y1 + 1;

        /* zeta = q'*x2 */
        zeta = 0;
        for (i=0; i<conesize-1; i++) {
            zeta += q[i]*x2[i];
        }

        /* y1 += eta^2*[ (a^2 + w)x1 + c*zeta ] */
        y1[0] += eta_square*( (a*a+w)*x1[0] + c*zeta );

        /* y2 += eta^2*[ (c*q*x1 + x2 + d*q*zeta ] */
        temp = c*x1[0] + d*zeta;
        for (i=0; i<conesize_m1; i++) {
            y2[i] += eta_square*( temp*q[i] + x2[i] );
        }

        /* prepare index for next cone */
        cone_start += conesize;
    }
#endif
#ifdef EXPCONE
    scaleToAddExpcone(y,x,C->expc,C->nexc,cone_start);
#endif
}

/**
 * Fast left-division by scaling matrix.
 * Returns z = W\lambda
 */
void unscale(pfloat* lambda, cone* C, pfloat* z)
{
	idxint i, j, l, cone_start;
	pfloat zeta, factor;

	/* LP cone */
	for( i=0; i < C->lpc->p; i++ ){ z[i] = SAFEDIV_POS(lambda[i], C->lpc->w[i]); }

	/* Second-order cone */
	cone_start = C->lpc->p;
	for( l=0; l < C->nsoc; l++ ){

		/* zeta = q'*lambda1 */
		zeta = 0;
		for( i=1; i < C->soc[l].p; i++ ){ zeta += C->soc[l].q[i-1] * lambda[cone_start + i]; }

		/* factor = -lambda0 + zeta / (1+a); */
		factor = -lambda[cone_start] + SAFEDIV_POS(zeta,(1+C->soc[l].a));

		/* second pass (on k): write out result */
		z[cone_start] = SAFEDIV_POS( (C->soc[l].a*lambda[cone_start] - zeta), C->soc[l].eta );
		for( i=1; i < C->soc[l].p; i++ ){
			j = cone_start+i;
			z[j] = SAFEDIV_POS( (lambda[j] + factor*C->soc[l].q[i-1]), C->soc[l].eta );
		}

		cone_start += C->soc[l].p;
	}
}



/**
 * Conic product, implements the "o" operator, w = u o v
 * and returns e'*w (where e is the conic 1-vector)
 */
pfloat conicProduct(pfloat* u, pfloat* v, cone* C, pfloat* w)
{
	idxint i, j, k, cone_start, conesize;
	pfloat u0, v0, mu;

    mu = 0;
	k=0;

	/* LP cone */
	for( i=0; i < C->lpc->p; i++ ){
        w[k] = u[i] * v[i];
        mu += w[k] < 0 ? -w[k] : w[k];
        k++;
    }

	/* Second-order cone */
	cone_start = C->lpc->p;
	for( i=0; i < C->nsoc; i++ ){
		conesize = C->soc[i].p;
		u0 = u[cone_start];
		v0 = v[cone_start];
		w[k] = eddot(conesize, u+cone_start, v+cone_start);
        mu += w[k] < 0 ? -w[k] : w[k];
        k++;
		for( j=1; j < conesize; j++ ){ w[k++] = u0*v[cone_start+j] + v0*u[cone_start+j]; }
		cone_start += conesize;
	}

    return mu;
}



/**
 * Conic division, implements the "\" operator, v = u \ w
 */
void conicDivision(pfloat* u, pfloat* w, cone* C, pfloat* v)
{
	idxint i, j, k, cone_start, conesize;
	pfloat rho, zeta, u0, w0, factor, temp;

	/* LP cone */
	for( i=0; i < C->lpc->p; i++ ){ v[i] = SAFEDIV_POS(w[i],u[i]); }

	/* Second-order cone */
	cone_start = C->lpc->p;
	for( i=0; i < C->nsoc; i++ ){
		conesize = C->soc[i].p;
		u0 = u[cone_start]; w0 = w[cone_start];
		rho = u0*u0;   zeta = 0;
		for( j=1; j < conesize; j++ ){
			k = cone_start+j;
			rho -= u[k]*u[k];
			zeta += u[k]*w[k];
		}
        temp = SAFEDIV_POS(zeta,u0) - w0;
        factor = SAFEDIV_POS(temp,rho);
        temp = u0*w0 - zeta;
		v[cone_start] = SAFEDIV_POS(temp,rho);
		for( j=1; j < conesize; j++ ){
			k = cone_start+j;
			v[cone_start+j] = factor*u[k] + SAFEDIV_POS(w[k],u0);
		}
		cone_start += C->soc[i].p;
	}
}

/*
 * Returns details on second order cone
 * Purpose: cleaner code
 */
void getSOCDetails(socone *soc, idxint *conesize, pfloat* eta_square, pfloat* d1, pfloat* u0, pfloat* u1, pfloat* v1, pfloat **q)
{
#if CONEMODE == 0
    *conesize = soc->p;
    *eta_square = soc->eta_square;
    *d1 = soc->d1;
    *u0 = soc->u0;
    *u1 = soc->u1;
    *v1 = soc->v1;
    *q = soc->q;
#endif
}


/*
 * Returns dx, dy and dz from the expanded and permuted version of
 * a search direction vector.
 */
void unstretch(idxint n, idxint p, cone *C, idxint *Pinv, pfloat *Px, pfloat *dx, pfloat *dy, pfloat *dz)
{
    idxint i,j,k,l;
    k = 0;
    for( i=0; i<n; i++ ){ dx[i] = Px[Pinv[k++]]; }
    for( i=0; i<p; i++ ){ dy[i] = Px[Pinv[k++]]; }
    j = 0;
    for( i=0; i<C->lpc->p; i++ ){ dz[j++] = Px[Pinv[k++]]; }
    for( l=0; l<C->nsoc; l++ ){
        for( i=0; i<C->soc[l].p; i++ ){ dz[j++] = Px[Pinv[k++]]; }
#if CONEMODE == 0
        k += 2;
#endif
    }
#ifdef EXPCONE
    for( l=0; l<C->nexc; l++)
    {
        for( i=0; i<3; i++ ){ dz[j++] = Px[Pinv[k++]]; }
    }
#endif
}
