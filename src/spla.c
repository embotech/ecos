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
 * Sparse matrix-transpose-vector multiplication.
 * 
 * 		y = X'*x
 *
 * where X is a sparse matrix and x and y are dense vectors.
 */
void sparseMtv(spmat* X, pfloat* x, pfloat* y){
	idxint i, j;
	for(j = 0; j < X->n; j++){
		/* Fill y with zeros */
		y[j] = 0;
		/* Add X'*x. Loop through columns instead of computing the transpose. */
		for(i = X->jc[j]; i < X->jc[j+1]; i++){
			y[j] += X->pr[i]*x[X->ir[i]];
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
 * Sparse matrix addition.
 * 
 *		S += X
 *
 * S needs to have enough space allocated.
 */
void sparseAdd(spmat* X, spmat* S){
	idxint i, j, k, l, p, h, t, s, nnz_in_column, nnz_new, n_moved, a;
	/* Iterate through columns */
	for(j = 0; j < X->n; j++){
		/* Count non-zeros in resulting column */
		nnz_in_column = (S->jc[j+1]) - (S->jc[j]) + (X->jc[j+1]) - (X->jc[j]); /* Assuming every entery is in a different row. Maximum number of resulting nonzeros */
		a = S->jc[j];
		for(i = X->jc[j]; i < X->jc[j+1]; i++){
			for(k = a; k < S->jc[j+1]; k++){
				/* If two entries are in the same place, reduce the number of nonzeros by one, and add the value of the entry in X to the entry in S */
				if(X->ir[i] == S->ir[k]){
					nnz_in_column--;
					a = k + 1;
					S->pr[k] += X->pr[i];
					break;
				}
			}
		}
		nnz_new = nnz_in_column - S->jc[j+1] + S->jc[j];
		if(nnz_new != 0){
		/* Move entries in ir and pr, change entries in jc */
			/* The entries in this and higher columns of S have to be moved to make space for the new entries from X*/ 
			for(l = S->jc[S->n]-1; l >= S->jc[j+1]-nnz_new; l--){
				if(l >= 0){
					S->ir[l+nnz_new] = S->ir[l];
					S->pr[l+nnz_new] = S->pr[l];
				}
			}
			/* The entries of the next columns now have a new index */
			for(l = j+1; l <= S->n; l++){
				/* New index: old one + new nonzeros before them */
				S->jc[l] += nnz_new;
			}
			/* Write new entries */
			a = S->jc[j];
			/* Count how many of the new values are yet written to S */
			n_moved = 0;
			/* Iterate through the current column */
			for(l = X->jc[j]; l < X->jc[j+1]; l++){
				if(n_moved == nnz_new){
					break;
				}
				t = X->ir[l]; /* The current row in X */
				/* Check if there are already entries in S in this column */
				if(S->jc[j+1]-nnz_new-S->jc[j] > 0){
					/* Iterate through the nonzeros in S and check where to put in the new entry coming from X */
					for(p = a; p < (S->jc[j+1])-nnz_new+n_moved; p++){
						h = S->ir[p]; /* The current row in S */
						/* If the entries are in the same row, they are already added up */
						if(t != h){
							if(p<S->jc[j+1]-nnz_new+n_moved-1){
								/* New entry comes in between values in the current column of S */
								if(t>h && t<S->ir[p+1]){
									for(s = S->jc[j+1]-nnz_new+n_moved; s > p+1; s--){
										S->ir[s] = S->ir[s-1];
										S->pr[s] = S->pr[s-1];
									}
									S->ir[p+1] = t;
									S->pr[p+1] = X->pr[l];
									n_moved++;
									
									break;
								}
							}
							/* New entry is the first entry in the column */
							else if(t<S->ir[a]){
								for(s = S->jc[j+1]-nnz_new+n_moved; s > a; s--){
									S->ir[s] = S->ir[s-1];
									S->pr[s] = S->pr[s-1];
								}
								S->ir[a] = t;
								S->pr[a] = X->pr[l];
								n_moved++;
								a++;								
								break;
							}
							/* New entry is the last entry in the column */
							else if(t>S->ir[S->jc[j+1]-nnz_new+n_moved-1]){
								S->ir[S->jc[j+1]-nnz_new+n_moved] = t;
								S->pr[S->jc[j+1]-nnz_new+n_moved] = X->pr[l];
								n_moved++;
								break;
							}
						}								
					}
				}
				/* If the current column of S has no entries, just write in the new entry */
				else{
					S->ir[S->jc[j]+n_moved] = t;
					S->pr[S->jc[j]+n_moved] = X->pr[l];
					n_moved++;
					
				}
			
			}
		}		
	}
}


/*
 * Sparse matrix-scalar division, S = S/eta.
 */
void sparseDivison(pfloat eta, spmat* S){
	idxint i;
	for(i = 0; i < S->nnz; i++){
		S->pr[i] = S->pr[i]/eta;
	}
}


/* 
 * Vector-scalar division x = x/eta. m = length(x).
 */
void vecDiv(pfloat eta, idxint m, pfloat* x){
	idxint i;
	if(eta > 0 && eta < EPS)
		eta = EPS;
	else if (eta < 0 && eta > -EPS)
		eta = -EPS;
	for(i = 0; i < m; i++){
		x[i] = x[i]/eta;
	}
}


/* 
 * Return product z = x'*y, where x and y are two columns of a sparse matrix, given as xir, xpr, xnnz, yir, ypr, ynnz 
 */
pfloat spmat_dotprod(idxint* xir, pfloat* xpr, idxint xnnz, idxint* yir, pfloat* ypr, idxint ynnz){
	pfloat z = 0;
	if(is_orthogonal(xir,yir,xnnz,ynnz)){
		return z;
	}
	idxint i, j;
	for(i = 0; i < xnnz; i++){
		for(j = 0; j < ynnz; j++){
			if(xir[i] == yir[j]){
				z += xpr[i]*ypr[j];
				break;
			}
		}
	}
	return z;	
}


/* 
 * Returns 1 if x and y are orthogonal, 0 otherwise.
 */
idxint is_orthogonal(idxint* x, idxint* y, idxint sizex, idxint sizey){
	idxint i, j;
	for(i = 0; i < sizex; i++){
		for(j = 0; j < sizey; j++){
			if(x[i]==y[j]){
				return 0;
			}
		}
	}
	return 1;
}


/*
 * For K >= 3 very accurate dot product x'*y. n; size of vectors. 
 */
pfloat DotK(idxint n, pfloat* x, pfloat* y, idxint K){
	idxint i;
	pfloat res, p, h, r[2*n], temp[2];
	TwoProduct(x[0],y[0],temp);
	p = temp[0];
	r[0] = temp[1];
	for(i = 1; i < n; i++){
		TwoProduct(x[i],y[i],temp);
		h = temp[0];
		r[i] = temp[1];
		TwoSum(p,h,temp);
		p = temp[0];
		r[n+i-1] = temp[1];
	}
	r[2*n-1] = p;
	res = SumK(2*n,r,K-1);
	return res;
}
/* Needed for DotK */
/* Computes the Vector sum. Higher K -> more accurate result */
pfloat SumK(idxint n, pfloat* p, idxint K){
	idxint k;
	pfloat res, temp[2];
	for(k = 0; k < K-1; k++){
		VecSum(n,p,temp);
	}
	res = 0;
	for(k = 0; k < n-1; k++){
		res += p[k];
	}
	res += p[n-1];
	return res;
}
/* Needed for SumK */
void VecSum(idxint n, pfloat* p, pfloat* vsum){
	idxint i;
	for(i = 1; i < n; i++){
		TwoSum(p[i],p[i-1],vsum);
		p[i] = vsum[0];
		p[i-1] = vsum[1];
	}
}
/* Accurate sum of two numbers. a+b = sum[0]+sum[1] */
void TwoSum(pfloat a, pfloat b, pfloat* sum){
	pfloat x, y, z;
	x = a+b;
	z = x-a;
	y = (a-(x-z))+(b-z);
	sum[0] = x;
	sum[1] = y;
}
/* Accurate product of two numbers. a*b = prod[0]+prod[1]. */ 
void TwoProduct(pfloat a, pfloat b, pfloat* prod){
	pfloat x, y, asplit[2], bsplit[2];
	x = a*b;
	Split(a,asplit);
	Split(b,bsplit);
	y = asplit[1]*bsplit[1]-(((x-asplit[0]*bsplit[0])-asplit[1]*bsplit[0])-asplit[0]*bsplit[1]);
	prod[0] = x;
	prod[1] = y;	
}
/* Needed for TwoProduct */
void Split(pfloat a, pfloat* split){
	pfloat c;
	pfloat factor = 134217729; /* = 2^s + 1, eps = 2^(-t), s = t/2, -> s = 27 in IEEE 754 double precision */
	c = factor*a;
	split[0] = c-(c-a);
	split[1] = a-split[0];
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
