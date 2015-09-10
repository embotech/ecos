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

/* Equilibration module (c) Eric Chu, March 2014 */

#include "ecos.h"

#if defined EQUILIBRATE && EQUILIBRATE > 0

#include <math.h>

#include "equil.h"
#include "spla.h"

void max_rows(pfloat *E, const spmat *mat)
{
    /* assumes mat is not null */
    idxint i, j, row;
    for(i = 0; i < mat->n; i++) { /* cols */
      for(j = mat->jc[i]; j < mat->jc[i+1]; j++) {
        row = mat->ir[j];
        E[row] = MAX(fabs(mat->pr[j]), E[row]);
      }
    }
}

void max_cols(pfloat *E, const spmat *mat)
{
    /* assumes mat is not null */
    idxint i, j;
    for(i = 0; i < mat->n; i++) { /* cols */
      for(j = mat->jc[i]; j < mat->jc[i+1]; j++) {
        E[i] = MAX(fabs(mat->pr[j]), E[i]);
      }
    }
}

void sum_sq_rows(pfloat *E, const spmat *mat)
{
    /* assumes mat is not null */
    idxint i, j, row;
    for(i = 0; i < mat->n; i++) { /* cols */
      for(j = mat->jc[i]; j < mat->jc[i+1]; j++) {
        row = mat->ir[j];
        E[row] += (mat->pr[j] * mat->pr[j]);
      }
    }
}

void sum_sq_cols(pfloat *E, const spmat *mat)
{
    /* assumes mat is not null */
    idxint i, j;
    for(i = 0; i < mat->n; i++) { /* cols */
      for(j = mat->jc[i]; j < mat->jc[i+1]; j++) {
        E[i] += (mat->pr[j] * mat->pr[j]);
      }
    }
}

void equilibrate_rows(const pfloat *E, spmat *mat)
{
    idxint i, j, row;

    for(i = 0; i < mat->n; i++) {
        /* equilibrate the rows of a matrix */
        for(j = mat->jc[i]; j < mat->jc[i+1]; j++) {
            row = mat->ir[j];
            mat->pr[j] /= E[row];
        }
    }
}

void equilibrate_cols(const pfloat *E, spmat *mat)
{
    idxint i, j;

    for(i = 0; i < mat->n; i++) {
        /* equilibrate the columns of a matrix */
        for(j = mat->jc[i]; j < mat->jc[i+1]; j++) {
            mat->pr[j] /= E[i];
        }
    }
}

void restore(const pfloat *D, const pfloat *E, spmat *mat)
{
    idxint i, j, row;

    for(i = 0; i < mat->n; i++) {
        /* equilibrate the rows of a matrix */
        for(j = mat->jc[i]; j < mat->jc[i+1]; j++) {
            row = mat->ir[j];
            mat->pr[j] *= (D[row] * E[i]);
        }
    }
}

void use_alternating_norm_equilibration(pwork *w)
{
    idxint i, j, ind;
    idxint num_cols = w->A ? w->A->n : w->G->n;
    idxint num_A_rows = w->A ? w->A->m : 0;
    idxint num_G_rows = w->G->m;
    pfloat sum;

    /* initialize equilibration vector to 0 */
    for(i = 0; i < num_cols; i++) {
        w->xequil[i] = 0.0;
    }
    for(i = 0; i < num_A_rows; i++) {
        w->Aequil[i] = 0.0;
    }
    for(i = 0; i < num_G_rows; i++) {
        w->Gequil[i] = 0.0;
    }

    /* compute norm across rows of A */
    if(w->A)
        sum_sq_rows(w->Aequil, w->A);
    /* compute norm across rows of G */
    if(num_G_rows > 0)
        sum_sq_rows(w->Gequil, w->G);

    /* now collapse cones together by taking average norm square */
    ind = w->C->lpc->p;
    for(i = 0; i < w->C->nsoc; i++) {
      sum = 0.0;
      for(j = 0; j < w->C->soc[i].p; j++) {
        sum += w->Gequil[ind + j];
      }
      for(j = 0; j < w->C->soc[i].p; j++) {
        w->Gequil[ind + j] = sum / w->C->soc[i].p;
      }
      ind += w->C->soc[i].p;
    }

#ifdef EXPCONE
    for(i = 0; i < w->C->nexc; i++) {
      sum = 0.0;
      for(j = 0; j < 3; j++) {
        sum += w->Gequil[ind + j];
      }
      for(j = 0; j < 3; j++) {
        w->Gequil[ind + j] = sum / 3.0;
      }
      ind += 3;
    }
#endif

    /* get the norm */
    for(i = 0; i < num_A_rows; i++) {
      w->Aequil[i] = fabs(w->Aequil[i]) < 1e-6 ? 1.0 : sqrt(w->Aequil[i]);
    }
    for(i = 0; i < num_G_rows; i++) {
      w->Gequil[i] = fabs(w->Gequil[i]) < 1e-6 ? 1.0 : sqrt(w->Gequil[i]);
    }

    /* now scale A */
    if(w->A)
        equilibrate_rows(w->Aequil, w->A);
    if(num_G_rows > 0)
        equilibrate_rows(w->Gequil, w->G);

    if(w->A)
        sum_sq_cols(w->xequil, w->A);
    if(num_G_rows > 0)
        sum_sq_cols(w->xequil, w->G);

    /* get the norm */
    for(i = 0; i < num_cols; i++) {
        w->xequil[i] = fabs(w->xequil[i]) < 1e-6 ? 1.0 : sqrt(w->xequil[i]);
    }
    if(w->A)
        equilibrate_cols(w->xequil, w->A);
    if(num_G_rows > 0)
        equilibrate_cols(w->xequil, w->G);

    /* the c vector is scaled in the ECOS_solve function
    for(i = 0; i < num_cols; i++) {
        w->c[i] /= w->xequil[i];
    }  */

    /* equilibrate the b vector */
    for(i = 0; i < num_A_rows; i++) {
        w->b[i] /= w->Aequil[i];
    }
    /* equilibrate the h vector */
    for(i = 0; i < num_G_rows; i++) {
        w->h[i] /= w->Gequil[i];
    }
}

void use_ruiz_equilibration(pwork *w)
{
    idxint i, j, ind, iter;
    idxint num_cols = w->A ? w->A->n : w->G->n;
    idxint num_A_rows = w->A ? w->A->m : 0;
    idxint num_G_rows = w->G->m;
    pfloat *xtmp = calloc(num_cols, sizeof(pfloat));
    pfloat *Atmp = calloc(num_A_rows, sizeof(pfloat));
    pfloat *Gtmp = calloc(num_G_rows, sizeof(pfloat));
    pfloat total;

    /* initialize equilibration vector to 1 */
    for(i = 0; i < num_cols; i++) {
        w->xequil[i] = 1.0;
    }
    for(i = 0; i < num_A_rows; i++) {
        w->Aequil[i] = 1.0;
    }
    for(i = 0; i < num_G_rows; i++) {
        w->Gequil[i] = 1.0;
    }

    /* iterative equilibration */
    for(iter = 0; iter < EQUIL_ITERS; iter++) {
        /* each iteration updates w->A and w->G */

        /* zero out the temp vectors */
        for(i = 0; i < num_cols; i++) {
            xtmp[i] = 0.0;
        }
        for(i = 0; i < num_A_rows; i++) {
            Atmp[i] = 0.0;
        }
        for(i = 0; i < num_G_rows; i++) {
            Gtmp[i] = 0.0;
        }

        /* compute norm across columns of A, G */
        if(w->A)
            max_cols(xtmp, w->A);
        if(num_G_rows > 0)
            max_cols(xtmp, w->G);

        /* compute norm across rows of A */
        if(w->A)
            max_rows(Atmp, w->A);

        /* compute norm across rows of G */
        if(num_G_rows > 0)
            max_rows(Gtmp, w->G);

        /* now collapse cones together by using total over the group */
        /* ECHU: not sure what the right thing to do here is */
        ind = w->C->lpc->p;
        for(i = 0; i < w->C->nsoc; i++) {
          total = 0.0;
          for(j = 0; j < w->C->soc[i].p; j++) {
            total += Gtmp[ind + j];
          }
          for(j = 0; j < w->C->soc[i].p; j++) {
            Gtmp[ind + j] = total;
          }
          ind += w->C->soc[i].p;
        }
#ifdef EXPCONE
       /*Do the same for the exponential cones*/
       for(i = 0; i < w->C->nexc; i++) {
         total = 0.0;
         for(j = 0; j < 3; j++) {
           total += Gtmp[ind + j];
         }
         for(j = 0; j < 3; j++) {
           Gtmp[ind + j] = total;
         }
         ind += 3;
       }
#endif


        /* take the sqrt */
        for(i = 0; i < num_cols; i++) {
          xtmp[i] = fabs(xtmp[i]) < 1e-6 ? 1.0 : sqrt(xtmp[i]);
        }
        for(i = 0; i < num_A_rows; i++) {
          Atmp[i] = fabs(Atmp[i]) < 1e-6 ? 1.0 : sqrt(Atmp[i]);
        }
        for(i = 0; i < num_G_rows; i++) {
          Gtmp[i] = fabs(Gtmp[i]) < 1e-6 ? 1.0 : sqrt(Gtmp[i]);
        }

        /* equilibrate the matrices */
        if(w->A)
            equilibrate_rows(Atmp, w->A);
        if(num_G_rows > 0)
            equilibrate_rows(Gtmp, w->G);

        if(w->A)
            equilibrate_cols(xtmp, w->A);
        if(num_G_rows > 0)
            equilibrate_cols(xtmp, w->G);

        /* update the equilibration matrix */
        for(i = 0; i < num_cols; i++) {
          w->xequil[i] *= xtmp[i];
        }
        for(i = 0; i < num_A_rows; i++) {
          w->Aequil[i] *= Atmp[i];
        }
        for(i = 0; i < num_G_rows; i++) {
          w->Gequil[i] *= Gtmp[i];
        }
    }

    /* the c vector is scaled in the ECOS_solve function
    for(i = 0; i < num_cols; i++) {
        w->c[i] /= w->xequil[i];
    } */

    /* equilibrate the b vector */
    for(i = 0; i < num_A_rows; i++) {
        w->b[i] /= w->Aequil[i];
    }
    /* equilibrate the h vector */
    for(i = 0; i < num_G_rows; i++) {
        w->h[i] /= w->Gequil[i];
    }

    free(xtmp);
    free(Atmp);
    free(Gtmp);
}

/* equilibrate */
void set_equilibration(pwork *w)
{
#if defined(RUIZ_EQUIL)
    use_ruiz_equilibration(w);
#elif defined(ALTERNATING_EQUIL)
    use_alternating_norm_equilibration(w);
#else
    /* use identity equilibration */
    idxint i;
    idxint num_cols = w->A ? w->A->n : w->G->n;
    idxint num_A_rows = w->A ? w->A->m : 0;
    idxint num_G_rows = w->G->m;

    /* initialize equilibration vector to 1 */
    for(i = 0; i < num_cols; i++) {
        w->xequil[i] = 1.0;
    }
    for(i = 0; i < num_A_rows; i++) {
        w->Aequil[i] = 1.0;
    }
    for(i = 0; i < num_G_rows; i++) {
        w->Gequil[i] = 1.0;
    }
#endif
}
/* invert the equilibration job */
void unset_equilibration(pwork *w)
{
    idxint i;
    /* idxint num_cols = w->A ? w->A->n : w->G->n; */
    idxint num_A_rows = w->A ? w->A->m : 0;
    idxint num_G_rows = w->G->m;

    if(w->A)
        restore(w->Aequil, w->xequil, w->A);
    if(num_G_rows > 0)
        restore(w->Gequil, w->xequil, w->G);

    /* the c vector is unequilibrated in the ECOS_solve function
        for(i = 0; i < num_cols; i++) {
        w->c[i] *= w->xequil[i];
    }*/

    /* unequilibrate the b vector */
    for(i = 0; i < num_A_rows; i++) {
        w->b[i] *= w->Aequil[i];
    }
    /* unequilibrate the h vector */
    for(i = 0; i < num_G_rows; i++) {
        w->h[i] *= w->Gequil[i];
    }
}

#endif
