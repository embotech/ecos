/*
 * TODO: MIT/BSD license this generated code
 *
 * This is ANSI C compatible template for the most generic code generation.
 * Only contains the parameter structure and the generated functions
 */

#ifndef __SQRT_H__
#define __SQRT_H__

#include "qcml_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the parameter struct */
typedef struct sqrt_params {

} sqrt_params;

/* the input dimensions struct */
typedef struct sqrt_dims {
    char SENTINEL; /* empty dims struct */
} sqrt_dims;

/* the solution struct
 *    users are responsible for keep track of the variable lengths
 */
typedef struct sqrt_sol {
    double x;
} sqrt_vars;
  
/* converts the 'sqrt' parameters into SOCP data
 *     allocates a qc_socp struct
 */
qc_socp * qc_sqrt2socp(const sqrt_params * params, const sqrt_dims * dims);

/* assigns the pointers for the variables in 'sqrt' to point to the proper
 * memory locations in the solution vector
 */
void qc_socp2sqrt(double * x, sqrt_vars * vars, const sqrt_dims * dims);

#ifdef __cplusplus
}
#endif
  
#endif // __SQRT_H__
