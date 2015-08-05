/**
 * Santiago Akle
 * ICME Stanford University 2014
 *
 * Computes the value \omega(z) defined as the solution y to
 * the equation y+log(y) = z for z real and z>=1.
 * Follows the recommendations by
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

#include "glblopts.h"

#if (defined _WIN32 || defined _WIN64 || defined _WINDLL )
#define _USE_MATH_DEFINES
#endif
#include <math.h>

pfloat wrightOmega(pfloat z);
