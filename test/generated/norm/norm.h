/*
 * TODO: MIT/BSD license this generated code
 *
 * This is ANSI C compatible template for the most generic code generation.
 * Only contains the parameter structure and the generated functions
 */

#ifndef __NORM_H__
#define __NORM_H__

#include "qcml_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the input dimensions struct */
typedef struct norm_dims {
    char SENTINEL; /* empty dims struct */
} norm_dims;

/* the solution struct
 *    users are responsible for keep track of the variable lengths
 */
typedef struct norm_vars {
    pfloat * x;
} norm_vars;
  
/* converts the 'norm' parameters into SOCP data
 *     allocates a qc_socp struct
 */
qc_socp * qc_norm2socp();

/* assigns the pointers for the variables in 'norm' to point to the proper
 * memory locations in the solution vector
 */
void qc_socp2norm(pfloat * x, norm_vars * vars, const norm_dims * dims);

#ifdef __cplusplus
}
#endif
  
#endif // __NORM_H__
