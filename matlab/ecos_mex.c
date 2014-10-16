/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2012-13 Alexander Domahidi [domahidi@control.ee.ethz.ch],
 * Automatic Control Laboratory, ETH Zurich.
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


#include "mex.h"
#include "matrix.h"
#include "ecos.h"
#include "ecos_bb.h"

/* THE mex-function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )  
{  
    idxint i, n, m, p;
    idxint exitcode;
    idxint numerr = 0;
    const mxArray* c;
    const mxArray* G;
    const mxArray* h;    
    const mxArray* A = NULL;    
    const mxArray* b = NULL;    
    const mxArray* dims;
    const mxArray* dims_l;
    const mxArray* dims_q;
    const mxArray* opts = NULL;
    const mxArray* opts_bool_idx = NULL;
    const mxArray* opts_verbose = NULL;
    const mxArray* opts_feastol = NULL;
    const mxArray* opts_reltol = NULL;
    const mxArray* opts_abstol = NULL;
    const mxArray* opts_feastol_inacc = NULL;
    const mxArray* opts_reltol_inacc = NULL;
    const mxArray* opts_abstol_inacc = NULL;
    const mxArray* opts_maxit = NULL;
    const mwSize *size_c;
    const mwSize *size_G;
    const mwSize *size_h;
    const mwSize *size_A;
    const mwSize *size_b;
    const mwSize *size_q;
    const mwSize* opts_bool_idx_size = NULL;
    
    const mwSize ZERO[2] = {0, 0};
    
    mxArray* outvar;
    mxArray* tinfos;

/* change number of infofields according to profiling setting */    
#if PROFILING > 0 
#define NINFOFIELDS 16
#else 
#define NINFOFIELDS 15
#endif
    const char *infofields[NINFOFIELDS] = { "exitflag",
                                            "infostring",
                                            "pcost",                                    
                                            "dcost",
                                            "pres",
                                            "dres",
                                            "pinf",
                                            "dinf",
                                            "pinfres",
                                            "dinfres",
                                            "gap",
                                            "relgap",   
                                            "r0",
                                            "numerr",
                                            "iter"                                      
#if PROFILING > 0                
                                           ,"timing"
#endif
                                           };
                                           
#if PROFILING == 1 
#define NTFIELDS 3    
    const char *tinfo[NTFIELDS] = {"runtime", "tsolve", "tsetup"};
#endif            
    
#if PROFILING == 2    
#define NTFIELDS 10
    const char *tinfo[NTFIELDS] = {"runtime", "tsolve", "tsetup", "tkktcreate", "tkktfactor", "tkktsolve", "torder", "ttranspose","tfactor_t1","tfactor_t2"};
#endif
                                   
            
    idxint l;
    double *q;
    idxint *qint;
    idxint *bool_vars_idx;
    idxint ncones;
    idxint numConicVariables = 0;
    idxint num_bool_vars;
    
    /* options */
    idxint verbose;
    pfloat abstol;
    pfloat reltol;
    pfloat feastol;
    idxint maxit;
    
    pfloat *Gpr = NULL;
    idxint *Gjc = NULL;
    idxint *Gir = NULL;
   
    pfloat *Apr = NULL;    
    idxint *Ajc = NULL;
    idxint *Air = NULL;   
    
    pfloat *cpr = NULL;
    pfloat *hpr = NULL;
    pfloat *bpr = NULL;
    
    pwork* mywork;
    ecos_bb_pwork* bb_pwork;
    

#ifdef MEXARGMUENTCHECKS     
    if( !(nrhs >= 4 && nrhs <= 7) )
    {
        mexPrintf("ECOS %s - (c) A. Domahidi, Automatic Control Laboratory, ETH Zurich, 2012-14.\n", ECOS_VERSION);
        mexErrMsgTxt("ECOS takes 4 to 7 arguments: ECOS(c,G,h,dims), ECOS(c,G,h,dims,opts), ECOS(c,G,h,dims,A,b), or ECOS(c,G,h,dims,A,b,opts)");
    }
#endif    
  
    /* get pointers to data */
    c = prhs[0];      size_c = c ? mxGetDimensions(c) : (const mwSize *) &ZERO;
    G = prhs[1];      size_G = G ? mxGetDimensions(G) : (const mwSize *) &ZERO;
    h = prhs[2];      size_h = h ? mxGetDimensions(h) : (const mwSize *) &ZERO;
    dims = prhs[3];    
    dims_l = dims ? mxGetField(dims, 0, "l") : NULL;
    dims_q = dims ? mxGetField(dims, 0, "q") : NULL; 
    size_q = dims_q ? mxGetDimensions(dims_q) : (const mwSize *) &ZERO;
    if( nrhs >= 6 )
    {
        A = prhs[4];  size_A = A ? mxGetDimensions(A) : (const mwSize *) &ZERO;
        b = prhs[5];  size_b = b ? mxGetDimensions(b) : (const mwSize *) &ZERO;
    }
    if( nrhs == 5 || nrhs == 7 )
    {
      opts = nrhs==5 ? prhs[4] : prhs[6]; 
      
      opts_bool_idx = opts ? mxGetField(opts, 0, "bool_vars_idx") : 0;
      if (opts_bool_idx != NULL){
        opts_bool_idx_size = mxGetDimensions(opts_bool_idx);
        /* Retrieve the number of boolean vars */
        num_bool_vars = (idxint) (opts_bool_idx_size[1] > opts_bool_idx_size[0] ? opts_bool_idx_size[1] : opts_bool_idx_size[0]);
      }
      
      opts_verbose = opts ? mxGetField(opts, 0, "verbose") : 0;
      opts_abstol = opts ? mxGetField(opts, 0, "abstol") : 0;
      opts_feastol = opts ? mxGetField(opts, 0, "feastol") : 0;
      opts_reltol = opts ? mxGetField(opts, 0, "reltol") : 0;
      opts_abstol_inacc = opts ? mxGetField(opts, 0, "abstol_inacc") : 0;
      opts_feastol_inacc = opts ? mxGetField(opts, 0, "feastol_inacc") : 0;
      opts_reltol_inacc = opts ? mxGetField(opts, 0, "reltol_inacc") : 0;
      opts_maxit = opts ? mxGetField(opts, 0, "maxit") : 0;
    }
    
    /* determine sizes */
    n = (idxint)size_c[0];
    m = (idxint)size_G[0];
    p = (nrhs >= 6) ? (idxint)size_A[0] : 0;
    
    /* argument checking */
#ifdef MEXARGMUENTCHECKS     
    if (opts_bool_idx != NULL){
        if (mxGetNumberOfDimensions(opts_bool_idx) != 2 ||
                (opts_bool_idx_size[1] != 1 && opts_bool_idx_size[0] != 1)){
            mexErrMsgTxt("Boolean index needs to be a vector");
        }        
    }
    
    if( !mxIsDouble(c) || size_c[1] > 1 )
    {
        mexErrMsgTxt("First argument needs to be a column vector");
    }

    if( !mxIsDouble(G) )
    {
        mexErrMsgTxt("Second argument needs to be a matrix");
    }    

    if( !mxIsSparse(G) )
    {
        mexErrMsgTxt("Second argument needs to be a sparse matrix");
    }    

    if( size_G[1] != n )
    {
        mexErrMsgTxt("c and G do not match in dimension");
    }

    if( m > 0 ){
        if( size_h[0] != m )
        {
            mexErrMsgTxt("G and h do not match in dimension");
        }
        
        if( size_h[1] != 1 )
        {
            mexErrMsgTxt("h is expected to be a column vector");
        }
    }

    if( !mxIsStruct(dims) )
    {
        mexErrMsgTxt("Struct dims expected as 4th argument");
    }

    if( dims_l == NULL && dims_q == NULL )	{
        mexErrMsgTxt("Neither dims.l nor dims.q exist - unconstrained problem?");
    } 
        
    if( mxIsSparse(c) )
    {
       mexErrMsgTxt("c must be a dense vector");
    } 

    if( mxIsSparse(h) )
    {
       mexErrMsgTxt("h must be a dense vector");
    } 

    if( nrhs >= 6 )
    {
        if( size_A[1] != n )
        {
            mexErrMsgTxt("c and A dimension mismatch");
        }
        
        if( !mxIsDouble(A) )
        {
            mexErrMsgTxt("Fifth argument needs to be a matrix");
        }    

        if( !mxIsSparse(A) )
        {
            mexErrMsgTxt("Fifth argument (A) needs to be a sparse matrix");
        } 
        
        if( size_b[0] != p )
        {
            mexErrMsgTxt("b and A dimension mismatch");
        }
        
        if( size_b[1] != 1 )
        {
            mexErrMsgTxt("b is expected to be a column vector");
        }
        
        if( mxIsSparse(b) )
        {
            mexErrMsgTxt("b must be a dense vector");
        } 
    }
    
    if( nrhs == 7 || nrhs == 5 )
    {
      if( !mxIsStruct(opts) )
      {
          mexErrMsgTxt("Struct opts expected as last argument");
      }
    }

    if( nlhs > 5 ){
         mexErrMsgTxt("ECOS has up to 5 output arguments only");
    }
#endif

    /* find out dimensions of cones */    
    l = (idxint)(*mxGetPr(dims_l)); numConicVariables += l;
    if( dims_q != NULL) {
        q = mxGetPr(dims_q);
    } else {
        q = NULL;
    }
    if( q == NULL || size_q[0] == 0 || size_q[1]==0 ){
        ncones = 0;
    } else {
        ncones = size_q[1] > size_q[0] ? size_q[1] : size_q[0];
    }
    
    /* get problem data in right format matrices */
    if( m > 0){
        Gpr = (pfloat *)mxGetPr(G);
        Gjc = (idxint *)mxGetJc(G);
        Gir = (idxint *)mxGetIr(G);
    } else { mexErrMsgTxt( "ECOS does not support equality constrained problems without inequalities yet." ); }
    if( p > 0 ){
        Apr = (pfloat *)mxGetPr(A);
        Ajc = (idxint *)mxGetJc(A);
        Air = (idxint *)mxGetIr(A);  
    }    
    cpr = (pfloat *)mxGetPr(c);
    hpr = (pfloat *)mxGetPr(h);
    if ( p > 0 ) {
      bpr = (pfloat *)mxGetPr(b);
    }
   
    
    /* we have to copy the number of cones since Matlab gives us double
     * values but the C version of ECOS expects idxint values */
    qint = (idxint *)mxMalloc(ncones*sizeof(idxint));
    for( i=0; i < ncones; i++ ){ qint[i] = (idxint)q[i]; numConicVariables += qint[i]; }
    
    
    /* check that number of rows matches info defined in dims */
    if ( numConicVariables != m ) {
        mexErrMsgTxt("Number of rows does not match information given in dims");
    }
    
    /* Switch between ecos_bb and ecos */
    if (opts_bool_idx != NULL){
        /* Shift the boolean indexes from matlab style (start at 1) to C style (start at 0) */
        bool_vars_idx = (idxint *)mxMalloc(num_bool_vars*sizeof(idxint));        
        for( i=0; i < num_bool_vars; ++i){ bool_vars_idx[i] = (idxint) ((mxGetPr(opts_bool_idx))[i] - 1); mexPrintf("\t%u : %u\n",i, bool_vars_idx[i]);}
        
        bb_pwork = ecos_bb_setup(n, m, p, l, ncones, qint, Gpr, Gjc, Gir, Apr, Ajc, Air, cpr, hpr, bpr, num_bool_vars, bool_vars_idx);
        
        mywork = bb_pwork->ecos_prob;

    }else{
        /* This calls ECOS setup function. */
        mywork = ECOS_setup(n, m, p, l, ncones, qint, Gpr, Gjc, Gir, Apr, Ajc, Air, cpr, hpr, bpr);
    }    
    
    if( mywork == NULL ){
        mexErrMsgTxt("Internal problem occurred in ECOS while setting up the problem.\nPlease send a bug report with data to Alexander Domahidi.\nEmail: domahidi@control.ee.ethz.ch");
    }
    
    /* Set options */
    if( opts != NULL ) {
        if(opts_verbose != NULL )
        {
            mywork->stgs->verbose = mxIsLogical(opts_verbose) ? (idxint) (*mxGetLogicals(opts_verbose)) : (idxint)(*mxGetPr(opts_verbose));
        }
        if(opts_feastol != NULL )
        {
            mywork->stgs->feastol = (pfloat)(*mxGetPr(opts_feastol));
        }
        if(opts_abstol != NULL )
        {
            mywork->stgs->abstol = (pfloat)(*mxGetPr(opts_abstol));
        }
        if(opts_reltol != NULL )
        {
            mywork->stgs->reltol = (pfloat)(*mxGetPr(opts_reltol));
        }
        if(opts_feastol_inacc != NULL )
        {
            mywork->stgs->feastol_inacc = (pfloat)(*mxGetPr(opts_feastol_inacc));
        }
        if(opts_abstol_inacc != NULL )
        {
            mywork->stgs->abstol_inacc = (pfloat)(*mxGetPr(opts_abstol_inacc));
        }
        if(opts_reltol_inacc != NULL )
        {
            mywork->stgs->reltol_inacc = (pfloat)(*mxGetPr(opts_reltol_inacc));
        }
        if(opts_maxit != NULL )
        {
            mywork->stgs->maxit = (idxint)(*mxGetPr(opts_maxit));
        }
        
    }
      
    /* Solve! */    
    /* Switch between ecos_bb and ecos */
    if (opts_bool_idx != NULL){
        exitcode = ecos_bb_solve(bb_pwork);
    }else{
        exitcode = ECOS_solve(mywork);
    }    
            
    /* create output */
    /* x */
    if( nlhs > 0 ){
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        mxSetPr(plhs[0], mywork->x);
        mxSetM(plhs[0], n);
        mxSetN(plhs[0], 1);
    }
        
    /* y */
    if( nlhs > 1 ){
        plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
        mxSetPr(plhs[1], mywork->y);
        mxSetM(plhs[1], mywork->p);
        mxSetN(plhs[1], 1);
    }
    
    if( nlhs > 2 ){
        plhs[2] = mxCreateStructMatrix(1, 1, NINFOFIELDS, infofields);
        
        /* 0. exitflag */
        outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)exitcode;
		mxSetField(plhs[2], 0, "exitflag", outvar);
        
        /* 1. primal objective */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->pcost;
		mxSetField(plhs[2], 0, "pcost", outvar);
        
        /* 2. dual objective */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->dcost;
		mxSetField(plhs[2], 0, "dcost", outvar);
        
        /* 3. primal residual */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->pres;
		mxSetField(plhs[2], 0, "pres", outvar);
        
        /* 4. dual residual */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->dres;
		mxSetField(plhs[2], 0, "dres", outvar);
        
        /* 5. primal infeasible? */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->pinf;
		mxSetField(plhs[2], 0, "pinf", outvar);
        
        /* 6. dual infeasible? */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->dinf;
		mxSetField(plhs[2], 0, "dinf", outvar);
        
        /* 7. primal infeasibility measure */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->pinfres;
		mxSetField(plhs[2], 0, "pinfres", outvar);
        
        /* 8. dual infeasibility measure */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->dinfres;
		mxSetField(plhs[2], 0, "dinfres", outvar);
        
         /* 9. duality gap */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->gap;
		mxSetField(plhs[2], 0, "gap", outvar);
        
        /* 10. relative duality gap */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->relgap;
		mxSetField(plhs[2], 0, "relgap", outvar);
        
        /* 11. feasibility tolerance??? */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->stgs->feastol;
		mxSetField(plhs[2], 0, "r0", outvar);
        
        /* 12. iterations */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->iter;
		mxSetField(plhs[2], 0, "iter", outvar);
        
        /* 13. infostring */
        if (opts_bool_idx != NULL){
            switch( exitcode ){
                case MI_OPTIMAL_SOLN:
                    outvar = mxCreateString("Optimal branch and bound solution found");
                    break;
                case MI_MAXITER_FEASIBLE_SOLN:
                    outvar = mxCreateString("Maximum iterations reached with feasible solution found");
                    break;
                case MI_MAXITER_NO_SOLN:
                    outvar = mxCreateString("Maximum iterations reached with no feasible solution found");
                    break;
                case MI_INFEASIBLE:
                    outvar = mxCreateString("Problem is infeasible");
                    break;
                default:
                    outvar = mxCreateString("UNKNOWN PROBLEM IN BRANCH AND BOUND SOLVER");
            }
        }else{
            switch( exitcode ){
                case ECOS_OPTIMAL:
                    outvar = mxCreateString("Optimal solution found");
                    break;
                case (ECOS_OPTIMAL + ECOS_INACC_OFFSET):
                    outvar = mxCreateString("Optimal solution found within reduced tolerances");
                    break;
                case ECOS_MAXIT:
                    outvar = mxCreateString("Maximum number of iterations reached");
                    break;
                case ECOS_PINF:
                    outvar = mxCreateString("Certificate of primal infeasibility found");
                    break;
                case (ECOS_PINF + ECOS_INACC_OFFSET):
                    outvar = mxCreateString("Certificate of primal infeasibility found within reduced tolerances");
                    break;
                case ECOS_DINF:
                    outvar = mxCreateString("Certificate of dual infeasibility found");
                    break;
                case (ECOS_DINF + ECOS_INACC_OFFSET):
                    outvar = mxCreateString("Certificate of dual infeasibility found within reduced tolerances");
                    break;
                case ECOS_NUMERICS:
                    outvar = mxCreateString("Numerical problems");
                    break;
                case ECOS_SIGINT:
                    outvar = mxCreateString("Interrupted by user (CTRL-C)");
                    break;
                case ECOS_OUTCONE:
                    outvar = mxCreateString("PROBLEM: Mulitpliers leaving the cone");
                    break;
                default:
                    outvar = mxCreateString("UNKNOWN PROBLEM IN SOLVER");
            }       
        }        
		mxSetField(plhs[2], 0, "infostring", outvar);
        
#if PROFILING > 0        
        /* 14. timing information */
		tinfos = mxCreateStructMatrix(1, 1, NTFIELDS, tinfo);
        
        /* 14.1 --> runtime */
        outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(outvar) = (double)mywork->info->tsolve + (double)mywork->info->tsetup;
		mxSetField(tinfos, 0, "runtime", outvar);
        
        /* 14.2 --> setup time */
        outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(outvar) = (double)mywork->info->tsetup;
		mxSetField(tinfos, 0, "tsetup", outvar);
        
        /* 14.3 --> solve time */
        outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(outvar) = (double)mywork->info->tsolve;
		mxSetField(tinfos, 0, "tsolve", outvar);

#if PROFILING > 1        
        
        /* 14.4 time to create KKT matrix */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->tkktcreate;
		mxSetField(tinfos, 0, "tkktcreate", outvar);

        /* 14.5 time for kkt solve */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->tkktsolve;
		mxSetField(tinfos, 0, "tkktsolve", outvar);
        
        /* 14.6 time for kkt factor */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->tfactor;
		mxSetField(tinfos, 0, "tkktfactor", outvar);
        
        /* 14.7 time for ordering */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->torder;
		mxSetField(tinfos, 0, "torder", outvar);
        
        /* 14.8 time for transposes */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->ttranspose;
		mxSetField(tinfos, 0, "ttranspose", outvar);
        
        /* 14.9 time for factoring, part 1 */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->tfactor_t1;
		mxSetField(tinfos, 0, "tfactor_t1", outvar);
        
        /* 14.10 time for transposes */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)mywork->info->tfactor_t2;
		mxSetField(tinfos, 0, "tfactor_t2", outvar);
#endif       
        
        mxSetField(plhs[2], 0, "timing", tinfos);        
#endif        
        
        /* 15. numerical error? */
        if( (exitcode == ECOS_NUMERICS) || (exitcode == ECOS_OUTCONE) || (exitcode == ECOS_FATAL) ){
            numerr = 1;
        }
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)numerr;
		mxSetField(plhs[2], 0, "numerr", outvar);        
    }
    
    /* s */
    if( nlhs > 3 ){
        plhs[3] = mxCreateDoubleMatrix(0, 0, mxREAL);
        mxSetPr(plhs[3], mywork->s);
        mxSetM(plhs[3], m);
        mxSetN(plhs[3], 1);
    }
    
    /* z */
    if( nlhs > 4 ){
        plhs[4] = mxCreateDoubleMatrix(0, 0, mxREAL);
        mxSetPr(plhs[4], mywork->z);
        mxSetM(plhs[4], m);
        mxSetN(plhs[4], 1);
    }
    
    /* cleanup */
    if (opts_bool_idx != NULL){
        ecos_bb_cleanup(bb_pwork, nlhs > 2? nlhs-1 : nlhs);
        mxFree(bool_vars_idx);
    }else{
        ECOS_cleanup(mywork, nlhs > 2? nlhs-1 : nlhs);    
    }
    
    mxFree(qint);
}