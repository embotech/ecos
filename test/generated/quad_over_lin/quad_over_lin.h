/*
 * TODO: MIT/BSD license this generated code
 *
 * This is ANSI C compatible template for the most generic code generation.
 * Only contains the parameter structure and the generated functions
 */

#ifndef __QUAD_OVER_LIN_H__
#define __QUAD_OVER_LIN_H__

#include "qcml_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the input dimensions struct */
typedef struct quad_over_lin_dims {
    char SENTINEL; /* empty dims struct */
} quad_over_lin_dims;

/* the solution struct
 *    users are responsible for keep track of the variable lengths
 */
typedef struct quad_over_lin_vars {
    pfloat * x;
} quad_over_lin_vars;
  
/* converts the 'quad_over_lin' parameters into SOCP data
 *     allocates a qc_socp struct
 */
qc_socp * qc_quad_over_lin2socp();

/* assigns the pointers for the variables in 'quad_over_lin' to point to the proper
 * memory locations in the solution vector
 */
void qc_socp2quad_over_lin(pfloat * x, quad_over_lin_vars * vars, const quad_over_lin_dims * dims);

#ifdef __cplusplus
}
#endif
  
#endif // __QUAD_OVER_LIN_H__
