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
 * The exponental cone module is (c) Santiago Akle, Stanford University,
 * [akle@stanford.edu]
 *
 */

#include "cone.h"
#include "glblopts.h"
#include "wright_omega.h"
#include <math.h>
#ifndef __ECOS_EXP_H__
#define __ECOS_EXP_H__

#if defined EXPCONE
/*
  * Exponential cone structure: We save the index where each expcone Hessian
  * column starts in the csr format.  For each cone we also allocate space for
  * the gradient of the barrier at each cone and the values of entries of the
  * Hessian matrix
  */
typedef struct expcone
{
    idxint colstart[3]; /* All cones are fixed size, we store the index
                         * where the column of the hessian starts in the
                         * permuted Newton system.
                         */
    pfloat v[6];        /* Uper triangular section of the hessian */
    pfloat g[3];        /* Gradient of the barrier */
} expcone;

/*
 * Evaluates the Hessian of the exponential dual cone barrier at the triplet
 * w[0],w[1],w[2], and stores the upper triangular part of the matrix mu*H(w)
 * at v[0],...,v[5]. The entries of the Hessian are arranged columnwise into v
 */
void evalExpHessian(pfloat* w, pfloat* v, pfloat mu);

/*
 * Evaluates the gradient of the dual exponential cone barrier g^\star(z) at the triplet
 * w[0],w[1],w[2], and stores the result at g[0],..,g[2].
 */
void evalExpGradient(pfloat* w, pfloat* g);

/*
 * Computes f_e(s_e) + f^\star_e(z_e)
 */
pfloat evalBarrierValue(pfloat* siter, pfloat *ziter, idxint fc, idxint nexc);

/*
 * Multiplies by y+=muH*x
 */
void scaleToAddExpcone(pfloat* y, pfloat* x, expcone* expcones, idxint nexc, idxint fc);

/*
 * Returns 1 if s is primal feasible w.r.t the primal exponential
 * cone and 0 i.o.c
 */
idxint evalExpPrimalFeas(pfloat *s, idxint nexc);

/*
 * Returns 1 if s is dual feasible w.r.t the dual exponential
 * cone and 0 i.o.c
 */
idxint evalExpDualFeas(pfloat *s, idxint nexc);

#endif
#endif /* End ifndef __ECOS_EXP_H__ */

