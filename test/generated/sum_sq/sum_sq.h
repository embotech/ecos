/*
 * TODO: MIT/BSD license this generated code
 *
 * This is ANSI C compatible template for the most generic code generation.
 * Only contains the parameter structure and the generated functions
 */

#ifndef __SUM_SQ_H__
#define __SUM_SQ_H__

#include "qcml_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the input dimensions struct */
typedef struct sum_sq_dims {
    char SENTINEL; /* empty dims struct */
} sum_sq_dims;

/* the solution struct
 *    users are responsible for keep track of the variable lengths
 */
typedef struct sum_sq_vars {
    pfloat * x;
} sum_sq_vars;
  
/* converts the 'sum_sq' parameters into SOCP data
 *     allocates a qc_socp struct
 */
qc_socp * qc_sum_sq2socp();

/* assigns the pointers for the variables in 'sum_sq' to point to the proper
 * memory locations in the solution vector
 */
void qc_socp2sum_sq(pfloat * x, sum_sq_vars * vars, const sum_sq_dims * dims);

#ifdef __cplusplus
}
#endif
  
#endif // __SUM_SQ_H__
