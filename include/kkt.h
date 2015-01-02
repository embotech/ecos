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


/* The KKT module.
 * Handles all computation related to KKT matrix:
 * - updating the matrix
 * - its factorization
 * - solving for search directions
 * - etc.
 */


#ifndef __KKT_H__
#define __KKT_H__

#include "glblopts.h"
#include "spla.h"
#include "cone.h"

typedef struct kkt{	
	spmat*  PKPt;    /* Permuted KKT matrix, upper part only      */	
	spmat*  L;       /* LDL factor L                              */
    
	pfloat* D;       /* diagonal matrix D                         */	
	pfloat* work1;   /* workspace needed for factorization        */
	pfloat* work2;   /* workspace needed for factorization        */
	pfloat* work3;   /* workspace needed for factorization        */
    pfloat* work4;   /* workspace needed for factorization        */
    pfloat* work5;   /* workspace needed for factorization        */
    pfloat* work6;   /* workspace needed for factorization        */
	pfloat* RHS1;    /* Right hand side 1						  */
	pfloat* RHS2;    /* Right hand side 2           			  */	
	pfloat* dx1;     /* search direction of size n				  */
	pfloat* dx2;     /* search direction of size n				  */
	pfloat* dy1;     /* search direction of size p				  */
	pfloat* dy2;     /* search direction of size p				  */
	pfloat* dz1;     /* search direction of size m				  */
	pfloat* dz2;     /* search direction of size m				  */	
	
    idxint* P;       /* permutation                               */
	idxint* Pinv;    /* reverse permutation						  */
	idxint* PK;      /* permutation of row indices of KKT matrix  */	
	idxint* Parent;  /* Elimination tree of factorization         */
	idxint* Sign;    /* Permuted sign vector for regularization   */
	idxint* Pattern; /* idxint workspace needed for factorization */
	idxint* Flag;    /* idxint workspace needed for factorization */
	idxint* Lnz;     /* idxint workspace needed for factorization */
	
	pfloat delta;    /* size of regularization					  */
} kkt;

/* Return codes */
#define KKT_PROBLEM (0)
#define KKT_OK      (1)

/* METHODS */

/** 
 * Factorization of KKT matrix. Just a convenient wrapper for the LDL call.
 * The second argument delta determindes the threshold of dynamic regularization,
 * while the last argument is the regularization parameter if it becomes active.
 *
 * If detailed profiling is turned on, the function returns the accumulated times
 * for sparsity pattern computation in t1 and for numerical solve in t2.
 */
#if PROFILING > 1
idxint kkt_factor(kkt* KKT, pfloat eps, pfloat delta, pfloat *t1, pfloat *t2);
#else
idxint kkt_factor(kkt* KKT, pfloat eps, pfloat delta);
#endif


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
idxint kkt_solve(kkt* KKT,
                 spmat* A, spmat* G,
                 pfloat* Pb,
                 pfloat* dx, pfloat* dy, pfloat* dz,
                 idxint n, idxint p, idxint m,
                 cone* C,
                 idxint isinit,
                 idxint nitref);


/**
 * Updates the permuted KKT matrix by copying in the new scalings.
 */
void kkt_update(spmat* PKP, idxint* P, cone *C);


/**
 * Initializes the (3,3) block of the KKT matrix to produce the matrix
 *
 * 		[0  A'  G']
 * K =  [A  0   0 ]
 *      [G  0  -I ]
 *
 * It is assumed that the A,G have been already copied in appropriately,
 * and that enough memory has been allocated (this is done in preproc.c module).
 *
 * Note that the function works on the permuted KKT matrix.
 */
void kkt_init(spmat* PKP, idxint* P, cone *C);


#endif
