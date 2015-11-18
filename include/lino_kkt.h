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

/********** FACTOR & SOLVE **********/

/* Factor */
void NEfactor(pfc* mypfc, cone* C);

/* Solve */
idxint NEsolve(pfc* mypfc,cone* C, pfloat* bx, pfloat* by, pfloat* bz);

#endif
