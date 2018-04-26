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
 * Interface for the built-in timer of ECOS.
 */
#ifndef __TIMER_H__
#define __TIMER_H__

#include "glblopts.h"

#if PROFILING > 0

#if (defined _WIN32 || defined _WIN64 || defined _WINDLL )

/* Use Windows QueryPerformanceCounter for timing */
#include <windows.h>

typedef struct timer{
	__int64 tic;
	__int64 toc;
	__int64 freq;
} timer;


#elif (defined __APPLE__)

#include <mach/mach_time.h>

/* Use MAC OSX  mach_time for timing */
typedef struct timer{
	uint64_t tic;
	uint64_t toc;
	mach_timebase_info_data_t tinfo;
} timer;



#else

/* Use POSIX clocl_gettime() for timing on non-Windows machines */
#include <time.h>
#include <sys/time.h>

typedef struct timer{
	struct timespec tic;
	struct timespec toc;
} timer;

#endif

/* METHODS are the same for both */
void tic(timer* t);
pfloat toc(timer* t);

#endif /* END IF PROFILING > 0 */

#endif
/* END IFDEF __TIMER_H__ */
