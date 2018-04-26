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
 * Interface for ECOS signal handling.
 *
 * This module is (c) Michael Grant, [mcg@cvxr.com] contributed by Github PR #82
 */

#ifndef __CTRLC_H__
#define __CTRLC_H__

#include "glblopts.h"

#if CTRLC > 0

/* METHODS are the same for both */
void init_ctrlc(void);
void remove_ctrlc(void);
int check_ctrlc(void);

#else /* CTRLC = 0 */

/* No signal handling. */
#define init_ctrlc()
#define remove_ctrlc()
#define check_ctrlc() 0

#endif /* END IF CTRLC > 0 */

#endif /* END IFDEF __TIMER_H__ */

