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


/* data type definitions used with ECOS */

#ifndef __GLBLOPTS_H__
#define __GLBLOPTS_H__

/* SET PRINT LEVEL ----------------------------------------------------- */
#define PRINTLEVEL (2)     /* 0: no prints					             */
						   /* 1: only final info				         */
                           /* 2: progress print per iteration            */
						   /* 3: debug level, enables print & dump fcns. */

#define MATLAB_FLUSH_PRINTS
                            /* print each iteration directly to Matlab.  */
                            /* this options considerably slows down the  */
                            /* solver, but is useful if you solve big    */
                            /* problems.                                 */

/* SET PROFILING LEVEL ------------------------------------------------- */
#define PROFILING (1)      /* 0: no timing information				     */
                           /* 1: runtime (divided in setup and solve)    */
                           /* 2: detailed profiling                      */

/* SET DEBUG LEVEL ----------------------------------------------------- */
#define DEBUG (0)          /* 0: no debugging information                */
                           /* 1: debug info & dump intermediate results  */
                           /* (flag used only for development)           */

/* DATA TYPES ---------------------------------------------------------- */
#include <float.h>
#include <math.h>
/* NOTE: Currently, pfloat MUST be double for ecos */
typedef double pfloat;              /* for numerical values  */
#define ECOS_INFINITY   (DBL_MAX + DBL_MAX)
#define ECOS_NAN        (ECOS_INFINITY - ECOS_INFINITY)

#if defined(_MSC_VER) && (_MSC_VER < 1900)
/* this will also return true if x is nan, but we don't check that anyway */
#define isinf(x) (!_finite(x))
#endif

/* Exponential cone */
#define EXPCONE      /*When defined the exponential cone solver code is enabled*/

/* SYSTEM INCLUDES FOR PRINTING ---------------------------------------- */
#if PRINTLEVEL > 0
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define PRINTTEXT mexPrintf
// #elif defined PYTHON
// #include <Python.h>
// #define PRINTTEXT PySys_WriteStdout
#else
#define PRINTTEXT printf
#endif
#include <stdio.h>
#else
#define PRINTTEXT(...)
#endif

#include "SuiteSparse_config.h"

/* USE SAME NUMBER REPRESENTATION FOR INDEXING AS AMD AND LDL ---------- */
#if defined(DLONG) && !defined(LDL_LONG) || defined(LDL_LONG) && !defined(DLONG)
#error "Inconsistent definition of DLONG and LDL_LONG"
#endif

#ifdef DLONG
typedef SuiteSparse_long idxint;
#else
typedef int idxint;
#endif

/* SYSTEM INCLUDE IF COMPILING FOR MATLAB ------------------------------ */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/* CHOOSE RIGHT MEMORY MANAGER ----------------------------------------- */
#ifdef MATLAB_MEX_FILE
#define MALLOC mxMalloc
#define FREE mxFree
#else
#define MALLOC malloc
#define FREE free
#endif

#endif
