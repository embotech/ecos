/*
 * TODO: MIT/BSD license this generated code
 *
 * This is ANSI C compatible template for the most generic code generation.
 * Only contains the parameter structure and the generated functions
 */

#ifndef __SQ_NORM_H__
#define __SQ_NORM_H__

#include "qcml_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the input dimensions struct */
typedef struct sq_norm_dims {
    char SENTINEL; /* empty dims struct */
} sq_norm_dims;

/* the solution struct
 *    users are responsible for keep track of the variable lengths
 */
typedef struct sq_norm_vars {
    pfloat * x;
} sq_norm_vars;
  
/* converts the 'sq_norm' parameters into SOCP data
 *     allocates a qc_socp struct
 */
qc_socp * qc_sq_norm2socp();

/* assigns the pointers for the variables in 'sq_norm' to point to the proper
 * memory locations in the solution vector
 */
void qc_socp2sq_norm(pfloat * x, sq_norm_vars * vars, const sq_norm_dims * dims);

#ifdef __cplusplus
}
#endif
  
#endif // __SQ_NORM_H__
