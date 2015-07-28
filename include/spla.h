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
 * Sparse linear algebra library for solver, i.e. no memory manager
 * such as malloc is accessed by this module.
 */

#ifndef __SPLA_H__
#define __SPLA_H__

#include "glblopts.h"

/* Data structure for sparse matrices */
typedef struct spmat{    
   idxint* jc;
   idxint* ir;
   pfloat* pr;
   idxint n;
   idxint m;
   idxint nnz;   
} spmat;


/* SPARSE MATRIX OPERATIONS PROVIDED BY THIS MODULE -------------------- */

/* 
 * Sparse matrix-vector multiply for operations
 * 
 *		y  =  A*x  (if a > 0 && newVector == 1)
 *		y +=  A*x  (if a > 0 && newVector == 0)
 *		y  = -A*x  (if a < 0 && newVector == 1)
 *		y -=  A*x  (if a < 0 && newVector == 0)
 *
 * where A is a sparse matrix and both x and y are assumed to be dense.
 */
void sparseMV(spmat* A, pfloat* x, pfloat* y, idxint a, idxint newVector);


/*
 * Sparse matrix-transpose-vector multiplication.
 * 
 * 		y = X'*x
 *
 * where X is a sparse matrix and x and y are dense vectors.
 */
void sparseMtv(spmat* X, pfloat* x, pfloat* y);


/* 
 * Sparse matrix-transpose-vector multiply with subtraction.
 *
 * If newVector > 0, then this computes y = -A'*x,
 *                           otherwise  y -= A'*x,
 *
 * where A is a sparse matrix and both x and y are assumed to be dense.
 * If skipDiagonal == 1, then the contributions of diagonal elements are
 * not counted.
 *
 * NOTE: The product is calculating without explicitly forming the
 *       transpose.
 */
void sparseMtVm(spmat* A, pfloat* x, pfloat* y, idxint newVector, idxint skipDiagonal);


/* 
 * Sparse matrix addition.
 * 
 *		S += X
 *
 * S needs to have enough space allocated.
 */
void sparseAdd(spmat* X, spmat* S);


/*
 * Sparse matrix-scalar division, S = S/eta.
 */
void sparseDivison(pfloat eta, spmat* S);


/* 
 * Vector-scalar division x = x/eta. m = length(x).
 */
void vecDiv(pfloat eta, idxint m, pfloat* x);


/* 
 * Return product z = x'*y, where x and y are two columns of a sparse matrix, given as xir, xpr, xnnz, yir, ypr, ynnz 
 */
pfloat spmat_dotprod(idxint* xir, pfloat* xpr, idxint xnnz, idxint* yir, pfloat* ypr, idxint ynnz);


/* 
 * Returns 1 if x and y are orthogonal, 0 otherwise.
 */
idxint is_orthogonal(idxint* x, idxint* y, idxint sizex, idxint sizey);


/*
 * For K >= 3 very accurate dot product x'*y. n; size of vectors. 
 */
pfloat DotK(idxint n, pfloat* x, pfloat* y, idxint K);


/*
 * Vector addition y += x of size n. 
 */
void vadd(idxint n, pfloat* x, pfloat* y);

/*
 * Vector subtraction with scaling: y -= a*x of size n. 
 */
void vsubscale(idxint n, pfloat a, pfloat* x, pfloat* y);

/*
 * 2-norm of a vector.
 */
pfloat norm2(pfloat* v, idxint n);

/*
 * inf-norm of a vector.
 */
pfloat norminf(pfloat* v, idxint n);

/*
 * ECOS dot product z = x'*y of size n.
 */
pfloat eddot(idxint n, pfloat* x, pfloat* y);


#endif 
