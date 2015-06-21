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
  * The exponential cone module is (c) Santiago Akle, Stanford University,
  * [akle@stanford.edu]
  */

#include "expcone.h"

#ifdef EXPCONE
/* Evaluates the Hessian of the exponential dual cone barrier at the triplet
 * w[0],w[1],w[2], and stores the upper triangular part of the matrix mu*H(w)
 * at v[0],...,v[5] where the entries are arranged columnwise
 */
void evalExpHessian(pfloat* w, pfloat* v, pfloat mu)
{
      /**
       *  l = log(-y/x);
       *  r = -x*l-x+w;
       *  He = [[ 1/x^2 - 1/(r*x) + l^2/r^2,           1/(r*y) + (l*x)/(r^2*y),     -l/r^2];
       *       [   1/(r*y) + (l*x)/(r^2*y), 1/y^2 - x/(r*y^2) + x^2/(r^2*y^2), -x/(r^2*y)];
       *       [                    -l/r^2,                        -x/(r^2*y),      1/r^2]];
       */

        pfloat x     = w[0];
        pfloat y     = w[1];
        pfloat z     = w[2];
        pfloat l     = log(-y/x);
        pfloat r     = -x*l-x+z;
        v[0]         = mu*((r*r-x*r+l*l*x*x)/(r*x*x*r));
        v[1]         = mu*((z-x)/(r*r*y));
        v[2]         = mu*((r*r-x*r+x*x)/(r*r*y*y));
        v[3]         = mu*(-l/(r*r));
        v[4]         = mu*(-x/(r*r*y));
        v[5]         = mu*(1/(r*r));
}

/* Evaluates the gradient of the exponential cone g(z) at the triplet
 * w[0],w[1],w[2], and stores the result at g[0],..,g[2]
 */
void evalExpGradient(pfloat* w, pfloat* g)
{
        pfloat x     = w[0];
        pfloat y     = w[1];
        pfloat z     = w[2];
        pfloat l     = log(-y/x);
        pfloat r     = -x*l-x+z;

        g[0]         = (l*x-r)/(r*x);
        g[1]         = (x-r)/(y*r);
        g[2]         = -1/r;
}

/* Computes f_e(s_e) + f^\star_e(z_e) */
pfloat evalBarrierValue(pfloat* siter, pfloat *ziter, idxint fc, idxint nexc)
{
    pfloat l, u, v, w, x, y, z, o;

    pfloat primal_barrier = 0.0;
    pfloat dual_barrier   = 0.0;

    idxint j;
    
	/* Move to the first exponential cone slack */
    ziter = ziter+fc;
    siter = siter+fc;
	
	/* For the dual cone measure -u,v, -ul-u+w */
    /* For the primal cone measure z,v,omega-1 */
    for(j=0;j<nexc;j++)
    {
        /* Extract the entries */
        u = ziter[0];
        v = ziter[1];
        w = ziter[2];

        x = siter[0];
        y = siter[1];
        z = siter[2];

        l = log(-v/u);
        dual_barrier += -log(w-u-u*l)-log(-u)-log(v);

        /* Primal Cone */
        o = wrightOmega(1-x/z-log(z)+log(y));
        o = (o-1)*(o-1)/o;
        primal_barrier += -log(o)-2*log(z)-log(y)-3;

        ziter += 3;
        siter += 3;

    }
    return primal_barrier+dual_barrier;
}


/*
 * Computes y[fc:end] += muH(x[fc:end])*x[fc:end], where
 * fc is the index of the first exponential slack.
 * This method assumes that the scalings have been updated by update scalings
 * and that C->expc[cone_number].v contains mu*H(x).
 *
 */
void scaleToAddExpcone(pfloat* y, pfloat* x, expcone* expc, idxint nexc, idxint fc)
{
    idxint l;
    /* Shift to the exponential slacks */
    x = x+fc;
    y = y+fc;

    for( l=0; l < nexc; l++ ){

        y[0]+= expc[l].v[0]*x[0]+expc[l].v[1]*x[1]+expc[l].v[3]*x[2];
        y[1]+= expc[l].v[1]*x[0]+expc[l].v[2]*x[1]+expc[l].v[4]*x[2];
        y[2]+= expc[l].v[3]*x[0]+expc[l].v[4]*x[1]+expc[l].v[5]*x[2];

        /* prepare index for next cone */
        x += 3;
        y += 3;
    }
}

/*
 * Returns 1 if s is primal feasible
 * with respect to the exponential cone,
 * and 0 i.o.c
 */
idxint evalExpPrimalFeas(pfloat *s, idxint nexc)
{
    pfloat x1,x2,x3,tmp1,psi;
    idxint j = 0;

    for(j =0 ; j < nexc; j++)
    {
       x1 = s[3*j];
       x2 = s[3*j+1];
       x3 = s[3*j+2];
       tmp1 = log(x2/x3);
       psi   = x3*tmp1 - x1;
       if(psi<0||x2<0||x3<0)
       {
            return 0;
       }

    }
    return 1;
}

/*
 * Returns 1 if s is dual feasible
 * with respect to the dual of the exponential cone,
 * and 0 i.o.c
 */
idxint evalExpDualFeas(pfloat *z, idxint nexc)
{
    pfloat x1,x2,x3,tmp1,psi;
    idxint j = 0;

    for(j =0 ; j < nexc; j++)
    {
       x1 = z[3*j];
       x2 = z[3*j+1];
       x3 = z[3*j+2];
       tmp1 = log(-x2/x1);
       psi   = -x1-x1*tmp1+x3;
       if(0<x1||x2<0||psi<0)
       {
            return 0;
       }

    }
    return 1;
}

#endif

