/* ========================================================================== */
/* === ldl.h:  include file for the LDL package ============================= */
/* ========================================================================== */

/* Copyright (c) Timothy A Davis, http://www.suitesparse.com.
 * All Rights Reserved.  See README for the License.
 *
 * Stripped down by Alexander Domahidi, 2012.
 */

#include "../../include/glblopts.h"
#include "../../include/ecos.h"

#include "SuiteSparse_config.h"

#ifdef LDL_LONG
#define LDL_int SuiteSparse_long
#define LDL_ID SuiteSparse_long_id

#define LDL_symbolic2 ldl_l_symbolic2
#define LDL_numeric2 ldl_l_numeric2
#define LDL_lsolve ldl_l_lsolve
#define LDL_lsolve2 ldl_l_lsolve2
#define LDL_dsolve ldl_l_dsolve
#define LDL_ltsolve ldl_l_ltsolve

#else
#define LDL_int int
#define LDL_ID "%d"

#define LDL_symbolic2 ldl_symbolic2
#define LDL_numeric2 ldl_numeric2
#define LDL_lsolve ldl_lsolve
#define LDL_lsolve2 ldl_lsolve2
#define LDL_dsolve ldl_dsolve
#define LDL_ltsolve ldl_ltsolve

#endif

/* ========================================================================== */
/* === int version ========================================================== */
/* ========================================================================== */

void ldl_symbolic2 (int n, int Ap [ ], int Ai [ ], int Lp [ ], int Parent [ ], int Lnz [ ], int Flag [ ]) ;

int ldl_numeric2 (int n, int Ap [ ], int Ai [ ], double Ax [ ],
    int Lp [ ], int Parent [ ], int Sign[], double eps, double delta, int Lnz [ ], int Li [ ], double Lx [ ],
    double D [ ], double Y [ ], int Pattern [ ], int Flag [ ]
#if PROFILING > 1
     ,double *t1, double *t2
#endif
                  ) ;

void ldl_lsolve (int n, double B [], int Lp [ ], int Li [ ], double Lx [ ]) ;
void ldl_lsolve2 (int n, double B [], int Lp [ ], int Li [ ], double Lx [ ], double X [ ]) ;

void ldl_dsolve (int n, double X [ ], double D [ ]) ;

void ldl_ltsolve (int n, double X [ ], int Lp [ ], int Li [ ],
    double Lx [ ]) ;

/* ========================================================================== */
/* === long version ========================================================= */
/* ========================================================================== */

void ldl_l_symbolic2 (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Lnz [ ],
    SuiteSparse_long Flag [ ]);

SuiteSparse_long ldl_l_numeric2 (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], double Ax [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Sign [ ], double eps, double delta, SuiteSparse_long Lnz [ ],
    SuiteSparse_long Li [ ], double Lx [ ], double D [ ], double Y [ ],
    SuiteSparse_long Pattern [ ], SuiteSparse_long Flag [ ]
#if PROFILING > 1
                                 ,double *t1, double *t2
#endif
                                 ) ;

void ldl_l_lsolve (SuiteSparse_long n, double B [ ], SuiteSparse_long Lp [ ], SuiteSparse_long Li [ ], double Lx [ ]) ;
void ldl_l_lsolve2 (SuiteSparse_long n, double B [ ], SuiteSparse_long Lp [ ], SuiteSparse_long Li [ ], double Lx [ ], double X [ ]) ;

void ldl_l_dsolve (SuiteSparse_long n, double X [ ], double D [ ]) ;

void ldl_l_ltsolve (SuiteSparse_long n, double X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], double Lx [ ]) ;

/* ========================================================================== */
/* === LDL version ========================================================== */
/* ========================================================================== */

#define LDL_DATE "April 6, 2013, with dynamic regularization by A. Domahidi"
#define LDL_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define LDL_MAIN_VERSION 2
#define LDL_SUB_VERSION 1
#define LDL_SUBSUB_VERSION 0
#define LDL_VERSION LDL_VERSION_CODE(LDL_MAIN_VERSION,LDL_SUB_VERSION)
