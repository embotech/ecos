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
	idxint p;         /* dimension of cone                               */
    pfloat eta;       /* eta                                             */
	pfloat etasqrt;   /* square root of eta                              */
    pfloat a;         /* a                                               */
	pfloat omega;     /* = q'*q                                          */
	pfloat atilde;    /* = a^2 + omega                                   */	
	//pfloat beta;      /* = (vtilde / qtilde)^2                         */
	idxint kkt_vU;    /* index for v in K (upper triangular part)        */	
	idxint kkt_qL;    /* index for q in K (lower triangular part)        */		
	idxint kkt_atilde;  /* index for atilde on diagonal                  */		
    idxint szidx;     /* index of associated slacks / multipliers        */
	pfloat* skbar;    /* temporary variables to work with                */
	pfloat* zkbar;    /* temporary variables to work with                */
	pfloat* q;        /* q vector (see documentation)                    */        
	pfloat* qtilde;   /* qtilde vector (see documentation)               */        
	pfloat* vtilde;   /* vtilde vector (see documentation)               */        
	idxint* kkt_qtU;  /* p-1 indices for q' in K (upper triangular part) */
    idxint* sidx;     /* p indices for sign entries in big sign vector   */
    pfloat qtildefact;
    pfloat vtildefact;
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

#endif