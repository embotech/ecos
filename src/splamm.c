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
 * Sparse linear algebra library for setup phase, i.e. this module
 * accesses MALLOC and hence should not go on an embedded platform.
 */

#include "splamm.h"

/* SYSTEM INCLUDES FOR MEMORY ALLOCATION ------------------------------- */
#if PRINTLEVEL > 2
#include <stdlib.h>
#endif


/* 
 * Cumulative sum. p[i] = sum w[0..i-1] and w = p on output.
 * Defined for vectors of length m.
 */
void spla_cumsum(idxint* p, idxint* w, idxint m)
{
	idxint cumsum = 0;
	idxint i;
	for( i=0; i < m; i++ ){ 
		p[i] = cumsum; 
		cumsum += w[i]; 
		w[i] = p[i];
	}
}


/**
 * Returns the inverse of permutation p of length n.
 */
void pinv(idxint n, idxint* p, idxint* pinv)
{
	idxint i;
	for( i=0; i<n; i++ ){ pinv[p[i]] = i; }
}


/**
 * Transpose a matrix; returns A = M',
 * and an index vector MtoMt that directly maps elements of M
 * to elements of M'.
 */
spmat* transposeSparseMatrix(spmat* M, idxint* MtoMt)
{	
	idxint j, i, k, q;    
	idxint* w;
  
  spmat* A = newSparseMatrix(M->n, M->m, M->nnz);
  if (M->nnz == 0) return A;
  
	w = (idxint *)MALLOC(M->m*sizeof(idxint));

	/* row count: how often does row k occur in M? */
	for( i=0; i < M->m; i++ ) { w[i] = 0; }
	for( k=0; k < M->nnz; k++ ) { w[M->ir[k]]++; }
	
	/* row pointers: cumulative sum of w gives A->jc */
	spla_cumsum(A->jc, w, M->m);

	/* now walk through M and copy data to right places and set row counter */
	for( j=0; j < M->n; j++ ){
		for( k = M->jc[j]; k < M->jc[j+1]; k++ ){
			q = w[M->ir[k]]++;
			A->ir[q] = j;
			A->pr[q] = M->pr[k];
			MtoMt[k] = q;
		}
	}

	FREE(w);
	return A;
}


/**
 * Create a new sparse matrix (uses MALLOC!) 
 */
spmat* newSparseMatrix(idxint m, idxint n, idxint nnz)
{
	idxint* jc = (idxint *)MALLOC((n+1)*sizeof(idxint));
	idxint* ir = (idxint *)MALLOC(nnz*sizeof(idxint));
	pfloat* pr = (pfloat *)MALLOC(nnz*sizeof(pfloat));
	jc[n] = nnz;
	return ecoscreateSparseMatrix(m, n, nnz, jc, ir, pr);
}


/** 
 * Create a sparse matrix from existing arrays (no MALLOC) 
 */
spmat* ecoscreateSparseMatrix(idxint m, idxint n, idxint nnz, idxint* jc, idxint* ir, pfloat* pr)
{
	spmat* M = (spmat *)MALLOC(sizeof(spmat));
	M->m = m;
	M->n = n;    
	M->nnz = nnz;	
	M->jc = jc;
    M->ir = ir;
    M->pr = pr;
	if (M->jc) M->jc[n] = nnz;
	return M;
}


/**
 * Create a new sparse matrix (uses FREE!)
 */
void freeSparseMatrix(spmat* M)
{
	if (M->ir) FREE(M->ir);
	if (M->jc) FREE(M->jc);
	if (M->pr) FREE(M->pr);
	FREE(M);
}


/**
 * Permutes a symmetric matrix with only the upper triangular part stored.
 * Writes the upper triangular part of C = A(p,p) in column compressed
 * storage format.
 *
 * The function additionally returns the mapping PK that maps the row indices
 * of the sparse matrix A on the row indices of C, such that C[P[k]] = A[k].
 *
 * NOTE: The matrix C and the vector PK are NOT created within this function 
 *       - you need to allocate them beforehand!!
 *
 * If PK is NULL then the last output argument is ignored.
 */
void permuteSparseSymmetricMatrix(spmat* A, idxint* pinv, spmat* C, idxint* PK)
{
	idxint i, i2, j, j2, k, q;
	idxint* w;

	/* create matrix, permutation vector and temporary work space */	
	w = (idxint *)MALLOC(A->n*sizeof(idxint));

	/* count number of entries per column of C, store this in w */
	for( j=0; j < A->n; j++ ){ w[j] = 0; }
	for( j=0; j < A->n; j++ ){
		j2 = pinv[j];
		for( k=A->jc[j]; k < A->jc[j+1]; k++ ){
			i = A->ir[k];
			if( i > j ) continue; /* skip lower triangular part of A */
			i2 = pinv[i];
			w[ i2 > j2 ? i2 : j2]++;
		}
	}

	/* cumulative sum gives column pointers */
	spla_cumsum(C->jc, w, A->n);

	/* copy data over */
	for( j=0; j < A->n; j++ ){
		j2 = pinv[j];
		for( k=A->jc[j]; k < A->jc[j+1]; k++ ){
			i = A->ir[k];
			if( i > j ) continue; /* skip lower triangular part of A */
			i2 = pinv[i];
			q = w[i2 > j2 ? i2 : j2]++;
			C->ir[q] = i2 < j2 ? i2 : j2;
			C->pr[q] = A->pr[k];
			if( PK ) PK[k] = q;
		}
	}

	/* free work space */
	FREE(w);
}


/*
 * Returns a copy of a sparse matrix A.
 */
spmat* copySparseMatrix(spmat* A)
{	
	idxint i;
	spmat* B = newSparseMatrix(A->m, A->n, A->nnz);

	/* copy over */
	for( i=0; i<=A->n; i++ ){ B->jc[i] = A->jc[i]; }
	for( i=0; i<A->nnz; i++ ){ B->ir[i] = A->ir[i]; }
	for( i=0; i<A->nnz; i++ ){ B->pr[i] = A->pr[i]; }

	return B;
}



/* ================================= DEBUG FUNCTIONS ======================= */
#if PRINTLEVEL > 2
/**
 * Prints a dense matrix.
 */
void printDenseMatrix(pfloat *M, idxint dim1, idxint dim2, char *name)
{
    idxint i,j;
    PRINTTEXT("%s = \n\t",name);
    for( i=0; i<dim1; i++ ){
        for( j=0; j<dim2; j++ ){
            if( j<dim2-1 )
                PRINTTEXT("% 14.12e,  ",M[i*dim2+j]);
            else
                PRINTTEXT("% 14.12e;  ",M[i*dim2+j]);
        }
        if( i<dim1-1){
            PRINTTEXT("\n\t");
        }
    }
    PRINTTEXT("\n");
}


/**
 * Prints a dense integer matrix.
 */
void printDenseMatrix_i(idxint *M, idxint dim1, idxint dim2, char *name)
{
    idxint i,j;
    PRINTTEXT("%s = \n\t",name);
    for( i=0; i<dim1; i++ ){
        for( j=0; j<dim2; j++ ){
            if( j<dim2-1 )
                PRINTTEXT("%d,  ",(int)M[i*dim2+j]);
            else
                PRINTTEXT("%d;  ",(int)M[i*dim2+j]);
        }
        if( i<dim1-1){
            PRINTTEXT("\n\t");
        }
    }
    PRINTTEXT("\n");
}


/*
 * Prints a sparse matrix.
 */
void printSparseMatrix(spmat* M)
{
    idxint j, i, row_strt, row_stop;
    idxint k = 0;
    for(j=0; j<M->n; j++){
        row_strt = M->jc[j];
        row_stop = M->jc[j+1];        
        if (row_strt == row_stop)
            continue;
        else {
            for(i=row_strt; i<row_stop; i++ ){                
                PRINTTEXT("\t(%3u,%3u) = %g\n", (int)M->ir[i]+1, (int)j+1, M->pr[k++]);
            }
        }
    }
}


/* dump a sparse matrix in Matlab format */
/* use LOAD and SPCONVERT to read in the file in MATLAB */
void dumpSparseMatrix(spmat* M, char* fn)
{	
	idxint j, i, row_strt, row_stop;
    idxint k = 0;
	FILE *f = fopen(fn,"w");
	if( f != NULL ){
		for(j=0; j<M->n; j++){
			row_strt = M->jc[j];
			row_stop = M->jc[j+1];        
			if (row_strt == row_stop)
				continue;
			else {
				for(i=row_strt; i<row_stop; i++ ){                
					fprintf(f,"%d\t%d\t%20.18e\n", (int)M->ir[i]+1, (int)j+1, M->pr[k++]);
				}
			}
		}
        fprintf(f,"%d\t%d\t%20.18e\n", (int)M->m, (int)M->n, 0.0);
		fclose(f);
		PRINTTEXT("File %s successfully written.\n", fn);
	} else {
		PRINTTEXT("Error during writing file %s.\n", fn);
	}
}


/**
 * Dumps a dense matrix of doubles to a CSV file.
 */
void dumpDenseMatrix(pfloat *M, int dim1, int dim2, char *fn)
{
    FILE *f = fopen(fn,"w");
    int i,j;
    
    /* swap dimensions if only vector to be printed */
    if( dim1 == 1 && dim2 > 1){
        i = dim2;
        dim2 = dim1;
        dim1 = i;
    }
    
    if( f != NULL ){
        for( i=0; i<dim1; i++ ){
            if( dim2 > 1 ){
                for( j=0; j<dim2-1; j++ ){
                    fprintf(f, "%+18.16e,",M[i*dim2+j]);
                }
                j = dim2;
                fprintf(f, "%+18.16e\n",M[i*dim2+j]);                
            } else {
                fprintf(f, "%+18.16e\n",M[i]);                
            }                
        }
        fclose(f);
        printf("Written %d x %d matrix to '%s'.\n",i,dim2,fn);
    } else {
        printf("ERROR: file %s could not be opened. Exiting.",fn);
        exit(1);
    }    
}


/**
 * Dumps a dense matrix of integers to a CSV file.
 */
void dumpDenseMatrix_i(idxint *M, int dim1, int dim2, char *fn)
{
    FILE *f = fopen(fn,"w");
    int i,j;
    
    /* swap dimensions if only vector to be printed */
    if( dim1 == 1 && dim2 > 1){
        i = dim2;
        dim2 = dim1;
        dim1 = i;
    }
    
    if( f != NULL ){
        for( i=0; i<dim1; i++ ){
            if( dim2 > 1 ){
                for( j=0; j<dim2-1; j++ ){
                    fprintf(f, "%d,",(int)M[i*dim2+j]);
                }
                j = dim2;
                fprintf(f, "%d\n",(int)M[i*dim2+j]);
            } else {
                fprintf(f, "%d\n",(int)M[i]);
            }                
        }
        fclose(f);
        printf("Written %d x %d matrix to '%s'.\n",i,dim2,fn);
    } else {
        printf("ERROR: file %s could not be opened. Exiting.",fn);
        exit(1);
    }    
}
#endif


