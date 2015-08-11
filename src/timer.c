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
 * Implements for the built-in timer of ECOS.
 *
 * Under Windows, we use QueryPerformanceCounter to obtain the counter value.
 *
 * Under Unix systems, we use clock_gettime function from <time.h>
 */


#include "timer.h"

#if PROFILING > 0

#if (defined WIN32 || _WIN64)

#include <windows.h>

void tic(timer* t)
{
	QueryPerformanceFrequency((LARGE_INTEGER*)&t->freq);
	QueryPerformanceCounter((LARGE_INTEGER*)&t->tic);
}

pfloat toc(timer* t)
{
	QueryPerformanceCounter((LARGE_INTEGER*)&t->toc);
	return ((t->toc - t->tic) / (pfloat)t->freq);
}


#elif (defined __APPLE__)

void tic(timer* t)
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

pfloat toc(timer* t)
{

    uint64_t duration; /* elapsed time in clock cycles*/
    
    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;
    
    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (pfloat)duration / 1000000000;
}



#else

/* read current time */
void tic(timer* t)
{
	clock_gettime(CLOCK_MONOTONIC, &t->tic);
}


/* return time passed since last call to tic on this timer */
double toc(timer* t)
{
	struct timespec temp;
    
	clock_gettime(CLOCK_MONOTONIC, &t->toc);	
    
	if ((t->toc.tv_nsec - t->tic.tv_nsec)<0) {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec-1;
		temp.tv_nsec = 1000000000+t->toc.tv_nsec - t->tic.tv_nsec;
	} else {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
		temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
	}
	return (pfloat)temp.tv_sec + (pfloat)temp.tv_nsec / 1000000000;
}

#endif


#endif
