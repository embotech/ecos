/*
 * Normal equations solver for ECOS */
 
#ifndef __lino_kkt_H__
#define __lino_kkt_H__

#include <stdio.h>
#include <math.h>

#include "cholmod.h"
#include "ecos.h"
#include "splamm.h"

/********** Product-form Cholesky struct **********/

typedef struct pfc{
	idxint ncones;
	pfloat delta;
	timer tnefactor, tnesolve;
	pfloat tfactor, tsolve;
	spmat *G, *A, *At;
	
	/* for product-form cholesky */
	spmat** G_br; /* blockrows of G */
	spmat** GtG; /* G_br'*G_br */
	spmat* S; /* sum(1/eta_square*GtG) */
	spmat* Spattern, Mpattern;
	pfloat** gtw; /* G_br'*wnew */
	pfloat** gte; /* G_br'*e0, e0 being the first unit vector */
	pfloat** wnew; /* new scaling */
	
	/* RHS stuff */
	pfloat* xpGtWinv2z;
	pfloat *bx, *by, *bz, *bztemp, *bxbybz;
	pfloat bxbybzsize;
	
	/* Solution */
	pfloat *dx, *dy, *dz, *workz;
	
	/* Iterative Refinement */
	pfloat *ex, *ey, *ez, *ddx, *ddy, *ddz;
	
	/* Cholmod */
	cholmod_common c;
	cholmod_sparse *Scm, *Scmreg, *Spatterncm, *Spatterncmreg, *Acm, *Atcm, *Gcm, *Z, *Zt, *M, *Mreg, *RegS, *RegM;
	cholmod_factor *L, *L_M;
	cholmod_dense *dxcm, *dycm, *workx, *worky, *xpGtWinv2zcm, *RHS, *RHStemp, *bzcm, *up_d, *down_d;
} pfc;

/********** INIT **********/

/* Setup, allocate memory. TESTED */
pfc* neSetup(idxint l, idxint ncones, idxint* q, spmat* G, spmat* A, pfloat delta);

/* Deallocate memory. TESTED */
void neCleanup(pfc* mypfc, idxint ncones, idxint l);

/* Compute amount of memory needed to store blockrow of X, row k+1 to (and including) row l+1. TESTED */
idxint mem_blockrow(spmat* X, idxint k, idxint l);

/* return blockrow of X, row k+1 to (and including) row l+1. uses MALLOC. TESTED */ 
spmat* blockrow(spmat* X, idxint k, idxint l);

/* Return Y = X'*X. uses MALLOC. TESTED */
spmat* sparseMtM(spmat* X);

/* Computes how much memory is needed to compute X'*X TESTED */
idxint count_mem(spmat* X);

/* Computes how much memory is needed to compute X'*X + eye*delta */
idxint count_mem_diag(spmat* X);

void initfactors(pfc* mypfc, cone* C);

/********** LinAlg **********/

/* write (k+1)-th column of sparse matrix X to x. TESTED */
void spmat_column(spmat* X, idxint k, pfloat* x);

/********** FACTOR & SOLVE **********/

/* Compute bxpGtWinv2bz (a part of the RHS) */
void xpGtWinv2z(pfc* mypfc, cone* C, idxint isItRef);

/* compute G'_i*w_i and G'_i*e0. TESTED */
void computeUpdates(pfc* mypfc, cone* C);

/* change scaling representation. TESTED */
void change_scaling(cone* C, pfloat** w_new);

/* Add the scaled G_i'*G_i matrices together */
void addS(pfc* mypfc, cone* C);

/* Regularize, then factor S = sum((1/eta^2_i)*G'_i*G_i) TESTED*/
void factorS(pfc* mypfc);

/* up- and downdates on factor L TESTED */
void updown(pfc* mypfc);

/* compute Z and M, L*Z = A', Z'*Z = M = A*(G'*W^(-2)*G)^(-1)*A' TESTED*/
void compZM(pfc* mypfc);

/* Regularize, then factor M = Z'*Z = A*(G'*W^(-2)*G)^(-1)*A' TESTED*/
void factorM(pfc* mypfc);

/* Compute RHS of normal equations form TESTED */
void RHS(pfc* mypfc, cone* C, idxint isItRef);

/* Solve TESTED */
void LinSyssolve(pfc* mypfc, cone* C, idxint isItRef);

/* Iterative refinement TESTED */
idxint itref(pfc* mypfc, cone* C);

/* Factor */
void NEfactor(pfc* mypfc, cone* C);

/* Solve */
idxint NEsolve(pfc* mypfc,cone* C, pfloat* bx, pfloat* by, pfloat* bz);

/********** DEBUG **********/

/* print sparse matrix */
void printSparse(spmat* Y);

/* print sparse cholmod matrix */
void printSparseCM(cholmod_sparse* Y, cholmod_common* c);


#endif
