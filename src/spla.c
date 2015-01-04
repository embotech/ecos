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

#include "spla.h"

#include <math.h>


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
void sparseMV(spmat* A, pfloat* x, pfloat* y, idxint a, idxint newVector)
{
	idxint i, j;
	/* fill y with zeros in case A happens to not have any zero rows */
	if( newVector > 0 ){ for( i=0; i<A->m; i++ ){ y[i] = 0; } }
	
  /* if there are no entries in A, just do nothing */
  if( A->nnz == 0 ) return;
  
	if( a > 0 ){
		/* Add A*x */
		for( j=0; j<A->n; j++ ){
			for( i=A->jc[j]; i < A->jc[j+1]; i++ ){
				y[A->ir[i]] += A->pr[i]*x[j];
			}
		}
	} else {
		/* Subtract A*x */
		for( j=0; j<A->n; j++ ){
			for( i=A->jc[j]; i < A->jc[j+1]; i++ ){
				y[A->ir[i]] -= A->pr[i]*x[j];
			}
		}
	}
}


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
void sparseMtVm(spmat* A, pfloat* x, pfloat* y, idxint newVector, idxint skipDiagonal)
{
	idxint i, j, k;
	/* fill y with zeros in case A happens to not have any zero rows */
	if( newVector > 0 ){ for( j=0; j<A->n; j++ ){ y[j] = 0; } }
  
  /* if there are no entries in A, just do nothing */
  if( A->nnz == 0 ) return;
	
	if( skipDiagonal ){
		/* Subtract A'*x */	
		for( j=0; j<A->n; j++ ){
			for( k=A->jc[j]; k < A->jc[j+1]; k++ ){		
				i = A->ir[k];
				y[j] -= i==j? 0 : A->pr[k]*x[i];
			}
		}
	} else {
		/* Subtract A'*x */	
		for( j=0; j<A->n; j++ ){
			for( k=A->jc[j]; k < A->jc[j+1]; k++ ){			
				y[j] -= A->pr[k]*x[A->ir[k]];
			}
		}
	}
}


/*
 * Vector addition y += x of size n. 
 */
void vadd(idxint n, pfloat* x, pfloat* y)
{
	idxint i;
	for( i=0; i<n; i++){ y[i] += x[i]; }
}


/*
 * Vector subtraction with scaling: y -= a*x of size n. 
 */
void vsubscale(idxint n, pfloat a, pfloat* x, pfloat* y)
{
	idxint i;
	for( i=0; i<n; i++){ y[i] -= a*x[i]; }
}

/*
 * 2-norm of a vector.
 */
pfloat norm2(pfloat* v, idxint n)
{
	idxint i;
	pfloat normsquare = 0;
	for( i=0; i<n; i++ ){ normsquare += v[i]*v[i]; }
	return sqrt(normsquare);
}


/*
 * infinity norm of a vector.
 */
pfloat norminf(pfloat* v, idxint n)
{
	idxint i;
	pfloat norm = 0;
    pfloat mv;
	for( i=0; i<n; i++ ){
        if( v[i] > norm ){ norm = v[i]; }
        mv = -v[i];
        if( mv > norm){ norm = mv; }
    }
	return norm;
}

/*
 * ECOS dot product z = x'*y of size n.
 */
pfloat eddot(idxint n, pfloat* x, pfloat* y)
{
	pfloat z = 0; idxint i;
	for( i=0; i<n; i++ ){ z += x[i]*y[i]; }
	return z;
}
