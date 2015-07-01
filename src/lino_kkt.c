#include "lino_kkt.h"
/* #include "mex.h" */

/********** INIT **********/

/* Setup, allocate memory */
pfc* neSetup(idxint l, idxint ncones, idxint* q, spmat* G, spmat* A, pfloat* bx, pfloat* by, pfloat* bz, pfloat delta){
	idxint i, conestart, S_nnz;
	/* Get pfc data structure */
	pfc* mypfc = (pfc*)MALLOC(sizeof(pfc));
	/* G */
	mypfc->G = G;
	/* A */
	mypfc->A = A;
	/* Allocate memory for pointers to matrices and vectors */
	S_nnz = count_mem_diag(G);
	mypfc->S = newSparseMatrix(G->n,G->n,S_nnz);
	for(i = 0; i <= G->n; i++){
		mypfc->S->jc[i] = 0;
	}
	for(i = 0; i < S_nnz; i++){
		mypfc->S->ir[i] = 0;
		mypfc->S->pr[i] = 0;
	}
	mypfc->GtG = (spmat**)MALLOC((l+ncones)*sizeof(spmat*));
	if(ncones){
		mypfc->G_br = (spmat**)MALLOC(ncones*sizeof(spmat*));
		mypfc->gtw = (pfloat**)MALLOC(ncones*sizeof(pfloat*));
		mypfc->gte = (pfloat**)MALLOC(ncones*sizeof(pfloat*));
		mypfc->wnew = (pfloat**)MALLOC(ncones*sizeof(pfloat*));
	}
	/* Allocate memory for blockrows&blockrow squares and compute them. Allocate memory for wnew, G'_i*w_i and G'_i*e0. */
	conestart = l; /* l = number of lp-cones */
	
	spmat* G_temp;
	for(i = 0; i < l; i++){
		G_temp = blockrow(G,i,i);
		mypfc->GtG[i] = sparseMtM(G_temp);
		freeSparseMatrix(G_temp);
	}
	for(i = 0; i < ncones; i++){	
		mypfc->G_br[i] = blockrow(G,conestart,conestart+q[i]-1);		
		mypfc->GtG[l+i] = sparseMtM(mypfc->G_br[i]);
		mypfc->gtw[i] = (pfloat*)MALLOC(G->n*sizeof(pfloat));
		mypfc->gte[i] = (pfloat*)MALLOC(G->n*sizeof(pfloat));
		mypfc->wnew[i] = (pfloat*)MALLOC(q[i]*sizeof(pfloat));
		conestart += q[i];
	}
	
	mypfc->ncones = ncones;
	mypfc->delta = delta;
	
	/* RHS */
	mypfc->xpGtWinv2z = (pfloat*)MALLOC(G->n*sizeof(pfloat));
	mypfc->bx = bx;
	mypfc->by = by;
	mypfc->bz = bz;
	mypfc->bxbybzsize = A->n+A->m+G->m;
	mypfc->bxbybz = (pfloat*)MALLOC(mypfc->bxbybzsize*sizeof(pfloat));
	
	/* Solution */
	mypfc->dx = (pfloat*)MALLOC(A->n*sizeof(pfloat));
	mypfc->dy = (pfloat*)MALLOC(A->m*sizeof(pfloat));
	mypfc->dz = (pfloat*)MALLOC(G->m*sizeof(pfloat));
	mypfc->workz = (pfloat*)MALLOC(G->m*sizeof(pfloat));
	
	/* Errors (for iterative refinement)*/
	mypfc->ex = (pfloat*)MALLOC(A->n*sizeof(pfloat));
	mypfc->ey = (pfloat*)MALLOC(A->m*sizeof(pfloat));
	mypfc->ez = (pfloat*)MALLOC((G->m)*sizeof(pfloat));
	
	/* Solution iterative refinement */
	mypfc->ddx = (pfloat*)MALLOC(A->n*sizeof(pfloat));
	mypfc->ddy = (pfloat*)MALLOC(A->m*sizeof(pfloat));
	mypfc->ddz = (pfloat*)MALLOC(G->m*sizeof(pfloat));	
	
	/* Cholmod stuff */
	cholmod_l_start(&(mypfc->c));
	
	mypfc->Scm = cholmod_l_allocate_sparse(mypfc->S->m,mypfc->S->n,mypfc->S->nnz,1,1,1,CHOLMOD_REAL,&(mypfc->c));		
	
	mypfc->Acm = cholmod_l_allocate_sparse(mypfc->A->m,mypfc->A->n,mypfc->A->nnz,1,1,0,CHOLMOD_REAL,&(mypfc->c));
	for(i = 0; i <= mypfc->A->n; i++){
		((idxint*)mypfc->Acm->p)[i] = mypfc->A->jc[i];
	}
	for(i = 0; i < mypfc->A->nnz; i++){
		((idxint*)mypfc->Acm->i)[i] = mypfc->A->ir[i];
		((pfloat*)mypfc->Acm->x)[i] = mypfc->A->pr[i];
	}
	
	mypfc->Atcm = cholmod_l_transpose(mypfc->Acm,2,&(mypfc->c));
	
	mypfc->Gcm = cholmod_l_allocate_sparse(mypfc->G->m,mypfc->G->n,mypfc->G->nnz,1,1,0,CHOLMOD_REAL,&(mypfc->c));
	for(i = 0; i <= mypfc->G->n; i++){
		((idxint*)mypfc->Gcm->p)[i] = mypfc->G->jc[i];
	}
	for(i = 0; i < mypfc->G->nnz; i++){
		((idxint*)mypfc->Gcm->i)[i] = mypfc->G->ir[i];
		((pfloat*)mypfc->Gcm->x)[i] = mypfc->G->pr[i];
	}
	
	mypfc->RegS = cholmod_l_speye(mypfc->G->n,mypfc->G->n,CHOLMOD_REAL,&(mypfc->c));
	mypfc->RegM = cholmod_l_speye(mypfc->A->m,mypfc->A->m,CHOLMOD_REAL,&(mypfc->c));
	
	mypfc->xpGtWinv2zcm = cholmod_l_allocate_dense(mypfc->A->n,1,mypfc->A->n,CHOLMOD_REAL,&(mypfc->c));	
	mypfc->RHS = cholmod_l_allocate_dense(mypfc->A->m,1,mypfc->A->m,CHOLMOD_REAL,&(mypfc->c));	
	mypfc->bzcm = cholmod_l_allocate_dense(mypfc->G->m,1,mypfc->G->m,CHOLMOD_REAL,&(mypfc->c));
	mypfc->up_d = cholmod_l_allocate_dense(mypfc->S->m,1,mypfc->S->m,CHOLMOD_REAL,&(mypfc->c));
	mypfc->down_d = cholmod_l_allocate_dense(mypfc->S->m,1,mypfc->S->m,CHOLMOD_REAL,&(mypfc->c));
	
	/* Cholmod settings */
	mypfc->c.final_ll = 0;
	mypfc->c.final_pack = 0;
	mypfc->c.supernodal = CHOLMOD_AUTO;
	/*
	mypfc->c.nmethods = 1;
	mypfc->c.method[0].ordering = CHOLMOD_NATURAL;	
	mypfc->c.postorder = 0;
	*/
	
	return mypfc;
}

/* Deallocate memory */
void neCleanup(pfc* mypfc, idxint ncones, idxint l){
	idxint i;
	/* for product-form cholesky */
	for(i = 0; i < l; i++){
		freeSparseMatrix(mypfc->GtG[i]);
	}	
	for(i = 0; i < ncones; i++){
		FREE(mypfc->gte[i]);
		FREE(mypfc->gtw[i]);
		freeSparseMatrix(mypfc->GtG[l+i]);
		freeSparseMatrix(mypfc->G_br[i]);
		FREE(mypfc->wnew[i]);
	}
	if(ncones){
		FREE(mypfc->gte);
		FREE(mypfc->gtw);
		FREE(mypfc->G_br);
		FREE(mypfc->wnew);
	}
	FREE(mypfc->GtG);
	freeSparseMatrix(mypfc->S);
	
	/* RHS */
	FREE(mypfc->xpGtWinv2z);
	FREE(mypfc->bxbybz);
	
	/* Solution */
	FREE(mypfc->dx);
	FREE(mypfc->dy);
	FREE(mypfc->dz);
	FREE(mypfc->workz);
	
	/* Iterative Refinement */
	FREE(mypfc->ex);
	FREE(mypfc->ey);
	FREE(mypfc->ez);
	FREE(mypfc->ddx);
	FREE(mypfc->ddy);
	FREE(mypfc->ddz);
	
	/* Cholmod stuff */ 
	cholmod_l_free_sparse(&(mypfc->Scm),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->Scmreg),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->Acm),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->Atcm),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->Gcm),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->RegS),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->RegM),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->M),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->Mreg),&(mypfc->c));
	cholmod_l_free_dense(&(mypfc->xpGtWinv2zcm),&(mypfc->c));
	cholmod_l_free_dense(&(mypfc->RHS),&(mypfc->c));
	cholmod_l_free_dense(&(mypfc->bzcm),&(mypfc->c));
	cholmod_l_free_dense(&(mypfc->up_d),&(mypfc->c));
	cholmod_l_free_dense(&(mypfc->down_d),&(mypfc->c));
	cholmod_l_free_factor(&(mypfc->L),&(mypfc->c));
	cholmod_l_free_factor(&(mypfc->L_M),&(mypfc->c));
	cholmod_l_finish(&(mypfc->c));
	
	FREE(mypfc);
}

/* Compute amount of memory needed to store blockrow of X, row k+1 to (and including) row l+1 */
idxint mem_blockrow(spmat* X, idxint k, idxint l){
	idxint i, j, p, mem_needed;
	mem_needed = 0;
	for(j = 0; j < X->n; j++){
		for(i = X->jc[j]; i < X->jc[j+1]; i++){
			for(p = k; p <= l; p++){
				if(X->ir[i] == p){
					mem_needed++;
					break;
				}
			}
		}
	}
	return mem_needed;
}

/* return blockrow of X, row k+1 to (and including) row l+1. Do this in init-phase! */
spmat* blockrow(spmat* X, idxint k, idxint l){
	idxint i, j, p, flag, nnz, mem_needed;
	mem_needed = mem_blockrow(X,k,l);
	spmat* B = newSparseMatrix(l-k+1,X->n,mem_needed);
	nnz = 0;
	for(j = 0; j < X->n; j++){
		flag = 0;
		for(i = X->jc[j]; i < X->jc[j+1]; i++){
			for(p = k; p <= l; p++){
				if(X->ir[i] == p){
					nnz += 1;
					if(flag == 0){
						B->jc[j] = nnz-1;
						flag = 1;
					}
					B->ir[nnz-1] = p-k;
					B->pr[nnz-1] = X->pr[i];	
				}
			}
		}
		if(!flag){
			B->jc[j] = nnz;
		}	
	}
	B->jc[B->n] = nnz;
	B->nnz = nnz;
	return B;
}

/* Return Y = X'*X */
spmat* sparseMtM(spmat* X){
	idxint i, j, k, l, mem, xnnz, ynnz, nnz, flag;
	nnz = 0;
	pfloat z;
	/* amount of memory needed? */
	mem = count_mem(X);
	/* allocate memory for result */
	spmat* Y = newSparseMatrix(X->n,X->n,mem);
	/* compute&fill in entries */
	for(j = 0; j < X->n; j++){
		Y->jc[j] = nnz;
		flag = 0;
		xnnz = (X->jc[j+1])-(X->jc[j]);
		if(xnnz != 0){
			idxint xir[xnnz];
			pfloat xpr[xnnz];
			for(i = X->jc[j]; i < X->jc[j+1]; i++){
				xir[i-(X->jc[j])] = X->ir[i];
				xpr[i-(X->jc[j])] = X->pr[i];	
			}
			for(k = 0; k < X->n; k++){
				ynnz = (X->jc[k+1]-X->jc[k]);
				if(ynnz != 0){
					idxint yir[ynnz];
					pfloat ypr[ynnz];
					for(l = X->jc[k]; l < X->jc[k+1]; l++){
						yir[l-(X->jc[k])] = X->ir[l];
						ypr[l-(X->jc[k])] = X->pr[l];
					}
					z = spmat_dotprod(xir,xpr,xnnz,yir,ypr,ynnz);
					if(z != 0){
						nnz += 1;
						if(!flag){
							Y->jc[j] = nnz-1;
							flag = 1;
						}
						Y->ir[nnz-1] = k;
						Y->pr[nnz-1] = z;
					}
				}
			}
		}
	}
	Y->jc[X->n] = nnz;
	Y->nnz = nnz;
	return Y;
}

/* Computes how much memory is needed to compute X'*X */
idxint count_mem(spmat* X){
	if(X->nnz == 0){
		return 0;
	}
	idxint i, j, k, l, sizex, sizey, mem_needed;
	/* amount of memory needed? */
	mem_needed = 0;
	for(j = 0; j < X->n; j++){
		sizex = (X->jc[j+1])-(X->jc[j]);
		if(sizex != 0){
			idxint x[sizex];
			for(i = X->jc[j]; i < X->jc[j+1]; i++){
				x[i-(X->jc[j])] = X->ir[i];
			}
			for(k = 0; k < X->n; k++){
				sizey = (X->jc[k+1])-(X->jc[k]);
				if(sizey != 0){
					idxint y[sizey];
					for(l = X->jc[k]; l < X->jc[k+1]; l++){
						y[l-(X->jc[k])] = X->ir[l];
					}
					if(!is_orthogonal(x,y,sizex,sizey)){
						mem_needed++;
					}
				}
			}
		}
	}
	return mem_needed;
}

/* Computes how much memory is needed to compute X'*X + eye*delta */
idxint count_mem_diag(spmat* X){
	if(X->nnz == 0){
		return 0;
	}
	idxint i, j, k, l, sizex, sizey, mem_needed, mem_diag;
	/* amount of memory needed? */
	mem_needed = 0;
	mem_diag = X->n;
	for(j = 0; j < X->n; j++){
		sizex = (X->jc[j+1])-(X->jc[j]);
		if(sizex != 0){
			idxint x[sizex];
			for(i = X->jc[j]; i < X->jc[j+1]; i++){
				x[i-(X->jc[j])] = X->ir[i];
			}
			for(k = 0; k < X->n; k++){
				sizey = (X->jc[k+1])-(X->jc[k]);
				if(sizey != 0){
					idxint y[sizey];
					for(l = X->jc[k]; l < X->jc[k+1]; l++){
						y[l-(X->jc[k])] = X->ir[l];
					}
					if(!is_orthogonal(x,y,sizex,sizey)){
						mem_needed++;
						if(j == k){
							mem_diag--;
						}
					}
				}
			}
		}
	}
	return mem_needed+mem_diag;
}


/********** LinAlg **********/

/* Sparse matrix addition, S += X. S needs to have enough space allocated. */
void sparseAdd(spmat* X, spmat* S){
	idxint i, j, k, l, p, h, t, s, nnz_in_column, nnz_new, n_moved, a;
	/* iterate through columns */
	for(j = 0; j < X->n; j++){
		/* count non-zeros in resulting column & add up entries in same row */
		nnz_in_column = (S->jc[j+1]) - (S->jc[j]) + (X->jc[j+1]) - (X->jc[j]);
		a = S->jc[j];
		for(i = X->jc[j]; i < X->jc[j+1]; i++){
			for(k = a; k < S->jc[j+1]; k++){
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
		/* move entries in ir and pr, change entries in jc */ 
			for(l = S->jc[S->n]-1; l >= S->jc[j+1]-nnz_new; l--){
				if(l >= 0){
					S->ir[l+nnz_new] = S->ir[l];
					S->pr[l+nnz_new] = S->pr[l];
				}
			}
			for(l = j+1; l <= S->n; l++){
				S->jc[l] += nnz_new;
			}
			/* write new entries */
			a = S->jc[j];
			n_moved = 0;
			for(l = X->jc[j]; l < X->jc[j+1]; l++){
				if(n_moved == nnz_new){
					break;
				}
				t = X->ir[l];
				if(S->jc[j+1]-nnz_new-S->jc[j] > 0){
					for(p = a; p < (S->jc[j+1])-nnz_new+n_moved; p++){
						h = S->ir[p];
						if(t != h){
							if(p<S->jc[j+1]-nnz_new+n_moved-1){
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
							else if(t>S->ir[S->jc[j+1]-nnz_new+n_moved-1]){
								S->ir[S->jc[j+1]-nnz_new+n_moved] = t;
								S->pr[S->jc[j+1]-nnz_new+n_moved] = X->pr[l];
								n_moved++;
								break;
							}
						}								
					}
				}
				else{
					S->ir[S->jc[j]+n_moved] = t;
					S->pr[S->jc[j]+n_moved] = X->pr[l];
					n_moved++;
					
				}
			
			}
		}		
	}
}

/* Sparse matrix-scalar division, S = S/eta. */
void sparseDivison(pfloat eta, spmat* S){
	idxint i;
	for(i = 0; i < S->nnz; i++){
		S->pr[i] = S->pr[i]/eta;
	}
}

/* Vector-scalar division x = x/eta. m = length(x). */
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

/* Sparse matrix-transpose-vector multiplication. y = X'*x */
void sparseMtv(spmat* X, pfloat* x, pfloat* y){
	idxint i, j;
	for(j = 0; j < X->n; j++){
		y[j] = 0;
		for(i = X->jc[j]; i < X->jc[j+1]; i++){
			y[j] += X->pr[i]*x[X->ir[i]];
		}
	}
}

/* Return product z = x'*y, where x and y are two columns of a sparse matrix, given as xir, xpr, xnnz, yir, ypr, ynnz */
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

/* write (k+1)-th column of sparse matrix X to x */
void spmat_column(spmat* X, idxint k, pfloat* x){
	idxint i;
	/* fill x with zeros */
	for(i = 0; i < X->m; i++){
		x[i] = 0;
	}
	for(i = X->jc[k]; i < X->jc[k+1]; i++){
		x[X->ir[i]] = X->pr[i];
	}
}

/* returns 1 if x and y are orthogonal, 0 otherwise. n-column vectors. */
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

/* Needed for TwoProduct */
void Split(pfloat a, pfloat* split){
	pfloat c;
	pfloat factor = 134217729; /* = 2^s + 1, eps = 2^(-t), s = t/2, -> s = 27 in IEEE 754 double precision */
	c = factor*a;
	split[0] = c-(c-a);
	split[1] = a-split[0];
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

/* Accurate sum of two numbers. a+b = sum[0]+sum[1] */
void TwoSum(pfloat a, pfloat b, pfloat* sum){
	pfloat x, y, z;
	x = a+b;
	z = x-a;
	y = (a-(x-z))+(b-z);
	sum[0] = x;
	sum[1] = y;
}

/* Needed for dotprod */
void FastTwoSum(pfloat a, pfloat b, pfloat* fsum){
	fsum[0] = a+b;
	fsum[1] = (a-fsum[0])+b;	
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

/* Accurate sum of a vector (1-norm) of length n */
pfloat Sum2s(idxint n, pfloat* p){
	idxint i;
	pfloat pi, sig, temp[2];
	pi = p[0];
	sig = 0;
	for(i = 1; i < n; i++){
		TwoSum(pi,p[i],temp);
		sig += temp[1];
		pi = temp[0];
	}
	return pi+sig;
}

/* Accurate dot-product x'*y */
pfloat Dot2s(idxint n, pfloat* x, pfloat* y){
	idxint i;
	pfloat p, s, h, r, q, temp[2];
	TwoProduct(x[0],y[0],temp);
	p = temp[0];
	s = temp[1];
	for(i = 1; i < n; i++){
		TwoProduct(x[i],y[i],temp);
		h = temp[0];
		r = temp[1];
		TwoSum(p,h,temp);
		p = temp[0];
		q = temp[1];
		s += (q+r);
	}
	return p+s;
}

/* Very accurate dot-product */
pfloat DotXBLAS(idxint n, pfloat* x, pfloat* y){
	idxint i;
	pfloat h, r, s, t, s1, t1, s2, t2, temp[2];
	s = 0;
	t = 0;
	for(i = 0; i < n; i++){
		TwoProduct(x[i],y[i],temp);
		h = temp[0];
		r = temp[1];
		TwoSum(s,h,temp);
		s1 = temp[0];
		s2 = temp[1];
		TwoSum(t,r,temp);
		t1 = temp[0];
		t2 = temp[1];
		s2 += t1;
		FastTwoSum(s1,s2,temp);
		t1 = temp[0];
		s2 = temp[1];
		t2 += s2;
		FastTwoSum(t1,t2,temp);
		s = temp[0];
		t = temp[1];
	}
	return s+t;
}

/* DotK: for K >= 3 very accurate dot product. n; size of vectors. */
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
/********** FACTOR & SOLVE **********/

/* Compute bxpGtWinv2bz (a part of the RHS) */
void xpGtWinv2z(pfc* mypfc, cone* C, idxint isItRef){
	idxint i;
	pfloat temp1[mypfc->G->m];
	pfloat temp2[mypfc->G->m];
	
	if(isItRef){
		unscale(mypfc->ez,C,temp1);
	}
	else{
		unscale(mypfc->bz,C,temp1);
	}
	unscale(temp1,C,temp2);
	sparseMtv(mypfc->G,temp2,mypfc->xpGtWinv2z);
	if(isItRef){
		vadd(mypfc->G->n,mypfc->ex,mypfc->xpGtWinv2z);
	}
	else{
		vadd(mypfc->G->n,mypfc->bx,mypfc->xpGtWinv2z);
	}
	/* Copy it to cholmod format */
	for(i = 0; i < mypfc->G->n; i++){
		((pfloat*)mypfc->xpGtWinv2zcm->x)[i] = mypfc->xpGtWinv2z[i];
	}
}

/* compute G'_i*w_i and G'_i*e0 */
void computeUpdates(pfc* mypfc, cone* C){
	idxint i, j, conestart;
	change_scaling(C,mypfc->wnew);
	conestart = 0;
	for(i = 0; i < C->nsoc; i++){
		pfloat e[C->soc[i].p];
		for(j = 0; j < C->soc[i].p; j++){
			e[j] = 0;
		}
		e[0] = SAFEDIV_POS(sqrt(2),C->soc[i].eta);
		vecDiv(C->soc[i].eta/sqrt(2),C->soc[i].p,mypfc->wnew[i]);
		sparseMtv(mypfc->G_br[i],mypfc->wnew[i],mypfc->gtw[i]);
		sparseMtv(mypfc->G_br[i],e,mypfc->gte[i]);
		conestart += C->soc[i].p;		
	}	
}

/* Change scaling representation for product form cholesky */
void change_scaling(cone* C,pfloat** w_new){
	idxint i, j, cone_start;
	cone_start = 0;
	/* SOC-cones (scaling stays the same for LP-cone) */
	for(i = 0; i < C->nsoc; i++){		
		pfloat a = C->soc[i].a; /* old wbar(1) */
		pfloat w = C->soc[i].w; /* old q'*q, where q = wbar(2:end) */ 
		pfloat* q = C->soc[i].q; /* old q = wbar(2:end) */ 
		idxint sizeq = (C->soc[i].p)-1; /* size of q */ 		
		pfloat aq[sizeq];
		pfloat qw[sizeq];		 
		
		/* compute new wbar(1) */
		w_new[i][0] = sqrt(0.5*(a*a+w+1));
		/* compute new wbar(2:end) */
		for(j = 0; j < sizeq; j++){
			aq[j] = -a*q[j]/(2*w_new[i][0]);
			qw[j] = w/(1+a)*q[j]+q[j];
		}
		vsubscale(sizeq,1/(2*w_new[i][0]),qw,aq);
		for(j = 0; j < C->soc[i].p-1; j++){
			w_new[i][j+1] = aq[j];
		}
		cone_start += C->soc[i].p;	
	}	
}

/* Add the scaled G_i'*G_i matrices together */
void addS(pfc* mypfc, cone* C){
	idxint i;
   	for(i = 0; i < C->lpc->p; i++){            
        sparseDivison(C->lpc->v[i],mypfc->GtG[i]);
        sparseAdd(mypfc->GtG[i],mypfc->S);
        
    }
        
    for(i = 0; i < C->nsoc; i++){
        sparseDivison(C->soc[i].eta_square,mypfc->GtG[i+C->lpc->p]);
        sparseAdd(mypfc->GtG[i+C->lpc->p],mypfc->S);
    }
}

/* Regularize, then factor S = sum((1/eta^2_i)*G'_i*G_i) */
void factorS(pfc* mypfc){
	idxint i;
	/* Copy S to cholmod format */	
	for(i = 0; i <= mypfc->S->n; i++){
		((idxint*)mypfc->Scm->p)[i] = mypfc->S->jc[i];
	}
	for(i = 0; i < mypfc->S->nnz; i++){
		((idxint*)mypfc->Scm->i)[i] = mypfc->S->ir[i];
		((pfloat*)mypfc->Scm->x)[i] = mypfc->S->pr[i];
	}

	pfloat alpha[2] = {1,0};
	pfloat beta[2] = {mypfc->delta,0};
	mypfc->Scmreg = cholmod_l_add(mypfc->Scm,mypfc->RegS,alpha,beta,1,1,&(mypfc->c));
	mypfc->Scmreg->stype = 1;
	
	mypfc->L = cholmod_l_analyze(mypfc->Scmreg,&(mypfc->c));
	cholmod_l_factorize(mypfc->Scmreg,mypfc->L,&(mypfc->c));
	/*
	idxint i;
	cholmod_sparse* Lsp = cholmod_l_factor_to_sparse(mypfc->L,&(mypfc->c));
	idxint col[Lsp->ncol];
	for(i = 0; i < Lsp->ncol; i++){
		col[i] = i;
	}	
	cholmod_sparse* Snew = cholmod_l_aat(Lsp,col,Lsp->ncol,1,&(mypfc->c));
	printSparseCM(Snew,&(mypfc->c));
	*/
}
	

/* up- and downdates on factor L */
void updown(pfc* mypfc){
	idxint i, j;
	cholmod_sparse* up;
	cholmod_sparse* down;
	/* up- and downdates */	
	for(i = 0; i < mypfc->ncones; i++){
		for(j = 0; j < mypfc->S->m; j++){
			/* Copy updates to cholmod format */
			((pfloat*)mypfc->up_d->x)[j] = mypfc->gtw[i][j];
			((pfloat*)mypfc->down_d->x)[j] = mypfc->gte[i][j];
		}
		up = cholmod_l_dense_to_sparse(mypfc->up_d,1,&(mypfc->c));
		down = cholmod_l_dense_to_sparse(mypfc->down_d,1,&(mypfc->c));
		cholmod_l_updown(1,up,mypfc->L,&(mypfc->c));
		cholmod_l_updown(0,down,mypfc->L,&(mypfc->c));
		
		cholmod_l_free_sparse(&up,&(mypfc->c));
		cholmod_l_free_sparse(&down,&(mypfc->c));
	}
	
	/* convert to LL' instead of LDL' */
	cholmod_l_change_factor(CHOLMOD_REAL,1,0,1,1,mypfc->L,&(mypfc->c));
	/*
	printf("IS_LL = %i\n",mypfc->L->is_ll);
	cholmod_sparse* Lsp = cholmod_l_factor_to_sparse(mypfc->L,&(mypfc->c));
	idxint col[Lsp->ncol];
	for(i = 0; i < Lsp->ncol; i++){
		col[i] = i;
	}	
	cholmod_sparse* Snew = cholmod_l_aat(Lsp,col,Lsp->ncol,1,&(mypfc->c));
	printSparseCM(Snew,&(mypfc->c));
	*/
}

/* compute Z and M, L*Z = A', Z'*Z = M = A*(G'*W^(-2)*G)^(-1)*A' */
void compZM(pfc* mypfc){
	/*
	idxint i;
	cholmod_sparse* Lsp = cholmod_l_factor_to_sparse(mypfc->L,&(mypfc->c));
	idxint col[Lsp->ncol];
	for(i = 0; i < Lsp->ncol; i++){
		col[i] = i;
	}	
	cholmod_sparse* Snew = cholmod_l_aat(Lsp,col,Lsp->ncol,1,&(mypfc->c));
	Snew->stype = 1;
	printf("Snew->stype = %i\n",Snew->stype);
	printSparseCM(Snew,&(mypfc->c));
	*/
	
	cholmod_sparse* Z = cholmod_l_spsolve(4,mypfc->L,mypfc->Atcm,&(mypfc->c));	
	cholmod_sparse* Zt = cholmod_l_transpose(Z,1,&(mypfc->c));
	mypfc->M = cholmod_l_ssmult(Zt,Z,1,1,1,&(mypfc->c));
	cholmod_l_free_sparse(&Z,&(mypfc->c));
	cholmod_l_free_sparse(&Zt,&(mypfc->c));
	/*
	cholmod_l_print_sparse(mypfc->M,"M",&(mypfc->c));
	printSparseCM(mypfc->M,&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->Z),&(mypfc->c));
	cholmod_l_free_sparse(&(mypfc->Zt),&(mypfc->c));
	*/
}

/* Regularize, then factor M = Z'*Z = A*(G'*W^(-2)*G)^(-1)*A' */
void factorM(pfc* mypfc){
	/* Regularize M */
	pfloat alpha[2] = {1,0};
	pfloat beta[2] = {mypfc->delta,0};
	mypfc->Mreg = cholmod_l_add(mypfc->M,mypfc->RegM,alpha,beta,1,1,&(mypfc->c));
	mypfc->Mreg->stype = 1;
	
	mypfc->L_M = cholmod_l_analyze(mypfc->Mreg,&(mypfc->c));
	cholmod_l_factorize(mypfc->Mreg,mypfc->L_M,&(mypfc->c));
	
	/*	
	cholmod_l_change_factor(CHOLMOD_REAL,1,0,1,1,mypfc->L_M,&(mypfc->c));
	
	idxint i;
	cholmod_sparse* LMsp = cholmod_l_factor_to_sparse(mypfc->L_M,&(mypfc->c));
	idxint col[LMsp->ncol];
	
	for(i = 0; i < LMsp->ncol; i++){
		col[i] = i;
	}
	cholmod_sparse* Mnew = cholmod_l_aat(LMsp,col,LMsp->ncol,1,&(mypfc->c));
	printSparseCM(Mnew,&(mypfc->c));
	*/
}

/* Compute RHS of normal equations form */
void RHS(pfc* mypfc, cone* C, idxint isItRef){
	idxint i;
	if(isItRef){
		for(i = 0; i < mypfc->A->m; i++){
			((pfloat*)mypfc->RHS->x)[i] = mypfc->ey[i];
		}	
	}
	else{
		for(i = 0; i < mypfc->A->m; i++){
			((pfloat*)mypfc->RHS->x)[i] = mypfc->by[i];
		}
	}
	/* Compute bxpGtWinv2bz */
	xpGtWinv2z(mypfc,C,isItRef);
	
	/* Compute RHS */
	mypfc->RHStemp = cholmod_l_solve(1,mypfc->L,mypfc->xpGtWinv2zcm,&(mypfc->c));	
	pfloat alpha[2] = {1,0};
	pfloat beta[2] = {-1,0};
	cholmod_l_sdmult(mypfc->Acm,0,alpha,beta,mypfc->RHStemp,mypfc->RHS,&(mypfc->c));
	cholmod_l_free_dense(&(mypfc->RHStemp),&(mypfc->c));
}

/* Solve */
void NEsolve(pfc* mypfc, cone* C, idxint isItRef){
	idxint i;
	
	/* Solve for dy */ 
	mypfc->worky = cholmod_l_solve(1,mypfc->L_M,mypfc->RHS,&(mypfc->c));
	
	/* Solve for dx */
	pfloat alpha[2] = {-1,0};
	pfloat beta[2] = {1,0};
	cholmod_l_sdmult(mypfc->Acm,1,alpha,beta,mypfc->worky,mypfc->xpGtWinv2zcm,&(mypfc->c));
	mypfc->workx = cholmod_l_solve(1,mypfc->L,mypfc->xpGtWinv2zcm,&(mypfc->c));

	/* Solve for dz */
	if(isItRef){
		for(i = 0; i < mypfc->G->m; i++){
			((pfloat*)mypfc->bzcm->x)[i] = mypfc->ez[i];
		}
	}
	else{
		for(i = 0; i < mypfc->G->m; i++){
			((pfloat*)mypfc->bzcm->x)[i] = mypfc->bz[i];
		}
	}		
	pfloat* bztemp = mypfc->bzcm->x;
	cholmod_l_sdmult(mypfc->Gcm,0,beta,alpha,mypfc->workx,mypfc->bzcm,&(mypfc->c));
	pfloat temp[mypfc->G->m];
	unscale(bztemp,C,temp);
	unscale(temp,C,mypfc->workz);
	
	if(isItRef){
		for(i = 0; i < mypfc->A->m; i++){
			mypfc->ddy[i] = ((pfloat*)mypfc->worky->x)[i];
		}
		for(i = 0; i < mypfc->A->n; i++){
			mypfc->ddx[i] =((pfloat*)mypfc->workx->x)[i];
		}
		for(i = 0; i < mypfc->G->m; i++){
			mypfc->ddz[i] = mypfc->workz[i];
		}
			
	}
	else{
		for(i = 0; i < mypfc->A->m; i++){
			mypfc->dy[i] = ((pfloat*)mypfc->worky->x)[i];
		}
		for(i = 0; i < mypfc->A->n; i++){
			mypfc->dx[i] = ((pfloat*)mypfc->workx->x)[i];
		}
		for(i = 0; i < mypfc->G->m; i++){
			mypfc->dz[i] = mypfc->workz[i];
		}		
	}
	cholmod_l_free_dense(&(mypfc->worky),&(mypfc->c));
	cholmod_l_free_dense(&(mypfc->workx),&(mypfc->c));
}

/* Iterative refinement */
void itref(pfc* mypfc, cone* C){
	idxint i, j;
	pfloat nex, ney, nez, nerr;
	pfloat nItref = 3;
	
	for(i = 0; i < mypfc->A->n; i++){
		mypfc->bxbybz[i] = mypfc->bx[i];
	}
	for(i = 0; i < mypfc->A->m; i++){
		mypfc->bxbybz[mypfc->A->n+i] = mypfc->by[i];
	}
	for(i = 0; i < mypfc->G->m; i++){
		mypfc->bxbybz[mypfc->A->n+mypfc->A->m+i] = mypfc->bz[i];
	}
	pfloat bnorm = 1+norminf(mypfc->bxbybz,mypfc->bxbybzsize);
	for(i = 0; i < nItref; i++){
		/* errors */
		
		/* ex */
		for(j = 0; j < mypfc->A->n; j++){
			mypfc->ex[j] = mypfc->bx[j];
		}
		sparseMtVm(mypfc->A,mypfc->dy,mypfc->ex,0,0);
		sparseMtVm(mypfc->G,mypfc->dz,mypfc->ex,0,0);
		
		/* ey */
		for(j = 0; j < mypfc->A->m; j++){
			mypfc->ey[j] = mypfc->by[j];
		}
		sparseMV(mypfc->A,mypfc->dx,mypfc->ey,-1,0);
		
		/* ez */
		for(j = 0; j < mypfc->G->m; j++){
			mypfc->ez[j] = mypfc->bz[j];
		}
		sparseMV(mypfc->G,mypfc->dx,mypfc->ez,-1,0);
		/* scale2addNE(mypfc->dz,mypfc->ez,C); */
		pfloat temp1[mypfc->G->m];
		pfloat temp2[mypfc->G->m];
		scale(mypfc->dz,C,temp1);
		scale(temp1,C,temp2);
		vadd(mypfc->G->m,temp2,mypfc->ez);
		
		
		/* DEBUG print errors
		for(j = 0; j < mypfc->A->n; j++){
			printf("ex[%i] = %f\n",j,mypfc->ex[j]);
		}
		for(j = 0; j < mypfc->A->m; j++){
			printf("ey[%i] = %f\n",j,mypfc->ey[j]);
		}
		for(j = 0; j < mypfc->G->m; j++){
			printf("ez[%i] = %f\n",j,mypfc->ez[j]);
		}
		*/
		
		/* maximum errors (infinity norm) */
		nex = norminf(mypfc->ex,mypfc->A->n);
		ney = norminf(mypfc->ey,mypfc->A->m);
		nez = norminf(mypfc->ez,mypfc->G->m);
		nerr = MAX(nex,ney);
		nerr = MAX(nerr,nez);
		
		/* continue? */
		if(nerr < LINSYSACC*bnorm){
			break;
		}
		
		/* RHS */
		RHS(mypfc,C,1);
		
		/* Solve new system */ 
		NEsolve(mypfc,C,1);
		
		/* Add to solution */
		vadd(mypfc->A->n,mypfc->ddx,mypfc->dx);
		vadd(mypfc->A->m,mypfc->ddy,mypfc->dy);
		vadd(mypfc->G->m,mypfc->ddz,mypfc->dz);
	}
	
}

/********** DEBUG **********/

/* print sparse matrix */
void printSparse(spmat* Y){
	idxint i;
	printf("jc: ");
	for(i = 0; i <= Y->n; i++){
		printf("%i ",Y->jc[i]);
	}
	printf("\nir: ");
	for(i = 0; i < Y->nnz; i++){
		printf("%i ",Y->ir[i]);
	}
	printf("\npr: ");
	for(i = 0; i < Y->nnz; i++){
		printf("%f ",Y->pr[i]);
	}
	printf("\n");
}

/* print sparse cholmod matrix */
void printSparseCM(cholmod_sparse* Y,cholmod_common* c){
	idxint Ynnz = cholmod_l_nnz(Y,c);
	idxint Ycol = Y->ncol;
	idxint* jc = Y->p;
	idxint* ir = Y->i;
	pfloat* pr = Y->x;
	
	idxint i;
	printf("jc: ");
	for(i = 0; i <= Ycol; i++){
		printf("%i ",jc[i]);
	}
	printf("\nir: ");
	for(i = 0; i < Ynnz; i++){
		printf("%i ",ir[i]);
	}
	printf("\npr: ");
	for(i = 0; i < Ynnz; i++){
		printf("%f ",pr[i]);
	}
	printf("\n");
}

int main(){
	pfloat x[5] = {0.123456789, -1.23456789, 3.456789, 0.0056789, 5.6789};
	pfloat* xptr = x;
	xptr++;
	pfloat res1 = Dot2s(4,xptr,xptr);
	pfloat res2 = eddot(4,xptr,xptr);
	printf("res1 = %f, res2 = %f\n",res1,res2);
}

