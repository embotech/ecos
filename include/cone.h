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


/* cone module */

#ifndef __CONE_H__
#define __CONE_H__

#include "glblopts.h"

#define CONEMODE (0) /* 0: expand to sparse cones (ECOS standard)       */
                     /* 1: dense cones (slow for big cones) */
                     /* 2: dense of fixed size */


/* LP CONE ------------------------------------------------------------- */
typedef struct lpcone{
    idxint p;         /* dimension of cone                               */    
    pfloat* w;        /* scalings                                        */
	pfloat* v;        /* = w^2 - saves p multiplications                 */	
	idxint* kkt_idx;  /* indices of KKT matrix to which scalings w^2 map */
} lpcone;


/* SECOND-ORDER CONE --------------------------------------------------- */
/* (all KKT indices are in compressed column format pointing into Kpr)   */
typedef struct socone{    
	idxint p;           /* dimension of cone                             */
	pfloat* skbar;      /* temporary variables to work with              */
	pfloat* zkbar;      /* temporary variables to work with              */
    pfloat a;           /* = wbar(1)                                     */
    pfloat d1;          /* first element of D                            */
    pfloat w;           /* = q'*q                                        */  
    pfloat eta;         /* eta = (sres / zres)^(1/4)                     */
    pfloat eta_square;  /* eta^2 = (sres / zres)^(1/2)                   */
    pfloat* q;          /* = wbar(2:end)                                 */
#if CONEMODE == 0
    idxint* Didx;       /* indices for D                                 */
    pfloat u0;          /* eta                                           */
    pfloat u1;          /* u = [u0; u1*q]                                */
    pfloat v1;          /* v = [0; v1*q]                                 */
#endif
#if CONEMODE > 0
    idxint* colstart;   /* colstart[n] gives index in KKT matrix where
                           the nth column of this scaling matrix in (3,3)
                           block starts                                  */
    pfloat c;           /* = 1 + a + w/(1+a)                             */
    pfloat d;           /* = 1 + 2/(1+a) + w/(1+a)^2                     */
#endif
    
} socone;


/* GENERAL STRUCTURE FOR A CONE ---------------------------------------- */
typedef struct cone{
	lpcone* lpc;	    /* LP cone					    */
	socone* soc;	    /* Second-Order cone            */
	idxint nsoc;        /* number of second-order cones */
} cone;


/* ERROR CODES --------------------------------------------------------- */
#define INSIDE_CONE (0)
#define OUTSIDE_CONE (1)


/* METHODS ------------------------------------------------------------- */
/**
 * Scales a conic variable such that it lies strictly in the cone.
 * If it is already in the cone, r is simply copied to s.
 * Otherwise s = r + (1+alpha)*e where alpha is the biggest residual.
 */
void bring2cone(cone* C, pfloat* r, pfloat* s);


/**
 * Update scalings.
 * Returns OUTSIDE_CONE as soon as any multiplier or slack leaves the cone,
 * as this indicates severe problems.
 */
idxint updateScalings(cone* C, pfloat* s, pfloat* z, pfloat* lambda);


/**
 * Fast multiplication by scaling matrix.
 * Returns lambda = W*z
 */
void scale(pfloat* z, cone* C, pfloat* lambda);


/**
 * Fast multiplication with V := W^2.
 * Computes y += W^2*x;
 */
void scale2add(pfloat *x, pfloat* y, cone* C);


/**
 * Fast left-division by scaling matrix.
 * Returns z = W\lambda
 */
void unscale(pfloat* lambda, cone* C, pfloat* z);


/**
 * Conic product, implements the "o" operator, w = u o v
 * and returns e'*w (where e is the conic 1-vector)
 */
pfloat conicProduct(pfloat* u, pfloat* v, cone* C, pfloat* w);


/**
 * Conic division, implements the "\" operator, w = u \ v
 */
void conicDivision(pfloat* u, pfloat* v, cone* C, pfloat* w);


/*
 * Returns details on second order cone
 * Purpose: cleaner code
 */
void getSOCDetails(socone *soc, idxint *conesize, pfloat* eta_square, pfloat* d1, pfloat* u0, pfloat* u1, pfloat* v1, pfloat **q);


/*
 * Returns dx, dy and dz from the expanded and permuted version of
 * a search direction vector.
 */
void unstretch(idxint n, idxint p, cone *C, idxint *Pinv, pfloat *Px, pfloat *dx, pfloat *dy, pfloat *dz);

#endif
