/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2012-14 Alexander Domahidi [domahidi@control.ee.ethz.ch],
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


/* ECOS TESTER MODULE */
/* THE CODE FOR MINIMAL UNIT TESTING HAS BEEN TAKEN FROM http://www.jera.com/techinfo/jtns/jtn002.html */
#include <stdio.h>

#include "minunit.h"
#include "ecos.h"

/* Include Tests */
#include "MPC/MPC01.h"
#include "MPC/MPC02.h"
#include "generated_tests/generated_tests.h"
#include "feas_prob/feas.h"

int tests_run = 0;

static char * all_tests() {
    mu_run_test(test_MPC01);
    mu_run_test(test_MPC02);
    mu_run_test(test_norm);
    mu_run_test(test_quad_over_lin);
    mu_run_test(test_sq_norm);
    mu_run_test(test_sum_sq);
    mu_run_test(test_inv_pos);
    /* mu_run_test(test_sqrt); */
    mu_run_test(test_feas);
    return 0;
}

int main(void) {
    char *result = all_tests();
    if (result != 0) {
        printf("%s\n", result);
    }
    else {
        printf("ALL TESTS PASSED\n");
    }
    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
