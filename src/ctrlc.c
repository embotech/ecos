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
 * Implements signal handling (ctrl-c) for ECOS.
 *
 * Under Windows, we use SetConsoleCtrlHandler.
 * Under Unix systems, we use sigaction.
 * For Mex files, we use utSetInterruptEnabled/utIsInterruptPending.
 *
 * This module is (c) Michael Grant, [mcg@cvxr.com] contributed by Github PR #82
 */

#include "ctrlc.h"

#if CTRLC > 0

#if defined MATLAB_MEX_FILE

 /* No header file available here; define the prototypes ourselves */
extern bool utIsInterruptPending(void);
extern bool utSetInterruptEnabled(bool);

static int istate;
void init_ctrlc(void)
{
    istate = (int)utSetInterruptEnabled(true);
}
void remove_ctrlc(void)
{
    utSetInterruptEnabled((bool)istate);
}
int check_ctrlc(void)
{
    return (int)utIsInterruptPending();
}

#elif defined _WIN32 || defined _WIN64

/* Use Windows SetConsoleCtrlHandler for signal handling */
#include <windows.h>

static int int_detected;
BOOL WINAPI handle_ctrlc(DWORD dwCtrlType)
{
    if (dwCtrlType != CTRL_C_EVENT) return FALSE;
    int_detected = 1;
    return TRUE;
}
void init_ctrlc(void)
{
    int_detected = 0;
    SetConsoleCtrlHandler( handle_ctrlc, TRUE );
}
void remove_ctrlc(void)
{
    SetConsoleCtrlHandler( handle_ctrlc, FALSE );
}
int check_ctrlc(void)
{
    return int_detected;
}

#else /* Unix */

/* Use POSIX clocl_gettime() for timing on non-Windows machines */

#include <signal.h>
static int int_detected;
struct sigaction oact;
void handle_ctrlc(int dummy)
{
    int_detected = dummy?dummy:-1;
}
void init_ctrlc(void)
{
    int_detected = 0;
    struct sigaction act;
    act.sa_flags = 0;
    sigemptyset(&act.sa_mask);
    act.sa_handler = handle_ctrlc;
    sigaction(SIGINT,&act,&oact);
}
void remove_ctrlc(void)
{
    struct sigaction act;
    sigaction(SIGINT,&oact,&act);
}
int check_ctrlc(void)
{
    return int_detected;
}

#endif /* END IF MATLAB_MEX_FILE / WIN32 */

#endif /* END IF CTRLC > 0 */

