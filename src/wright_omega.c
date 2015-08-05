#include "wright_omega.h"

/**
 * Santiago Akle
 * ICME Stanford University 2014
 *
 * Computes the value \omega(z) defined as the solution y to
 * the equation y+log(y) = z ONLY FOR z real and z>=1.
 *
 * Loosely follows the recommendations by
 * PIERS W. LAWRENCE, ROBERT M. CORLESS, and DAVID J. JEFFREY.
 * Published in:
 * Algorithm 917: Complex Double-Precision Evaluation of the Wright \omega Function
 * ACM Transactions on Mathematical Software (TOMS) TOMS Homepage table of contents archive
 * Volume 38 Issue 3, April 2012
 * Article No.	 20
 * Publication Date	2012-04-01 (yyyy-mm-dd)
 * Publisher	ACM New York, NY, USA
 * ISSN: 0098-3500 EISSN: 1557-7295 doi>10.1145/2168773.2168779
 */
pfloat wrightOmega(pfloat z)
{
    pfloat w  = 0.0;
    pfloat r  = 0.0;
    pfloat q  = 0.0;
    pfloat zi = 0.0;
    
	if(z<0.0) return -1; /* Fail if the input is not supported */
	
	if(z<1.0+M_PI)      /* If z is between 0 and 1+pi */
    {
        q = z-1;
        r = q;
        w = 1+0.5*r;
        r *= q;
        w += 1/16.0*r;
        r *= q;
        w -= 1/192.0*r;
        r *= q;
        w -= 1/3072.0*q;
        r *= q; /* (z-1)^5 */
        w += 13/61440.0*q;
        /* Initialize with the taylor series */
    }
    else
    {
        r = log(z);
        q = r;
        zi  = 1.0/z;
        w = z-r;
        q = r*zi;
        w += q;
        q = q*zi;
        w += q*(0.5*r-1);
        q = q*zi;
        w += q*(1/3.0*r*r-3.0/2.0*r+1);
        /* Initialize with w(z) = z-r+r/z^2(r/2-1)+r/z^3(1/3ln^2z-3/2r+1) */
    }
    /* FSC iteration */
    /* Initialize the residual */
    r = z-w-log(w);
    z = (1+w);
    q = z+2/3.0*r;
    w *= 1+r/z*(z*q-0.5*r)/(z*q-r);
    r = (2*w*w-8*w-1)/(72.0*(z*z*z*z*z*z))*r*r*r*r;
    /* Check residual */
    /*  if(r<1.e-16) return w; */ /*XXX Just do two rounds */
    z = (1+w);
    q = z+2/3.0*r;
    w *= 1+r/z*(z*q-0.5*r)/(z*q-r);
    r = (2*w*w-8*w-1)/(72.0*(z*z*z*z*z*z))*r*r*r*r;

    return w;
}
