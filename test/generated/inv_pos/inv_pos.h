/*
 * TODO: MIT/BSD license this generated code
 *
 * This is ANSI C compatible template for the most generic code generation.
 * Only contains the parameter structure and the generated functions
 */

#ifndef __INV_POS_H__
#define __INV_POS_H__

#include "qcml_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the input dimensions struct */
typedef struct inv_pos_dims {
    char SENTINEL; /* empty dims struct */
} inv_pos_dims;

/* the solution struct
 *    users are responsible for keep track of the variable lengths
 */
typedef struct inv_pos_vars {
    pfloat x;
} inv_pos_vars;
  
/* converts the 'inv_pos' parameters into SOCP data
 *     allocates a qc_socp struct
 */
qc_socp * qc_inv_pos2socp();

/* assigns the pointers for the variables in 'inv_pos' to point to the proper
 * memory locations in the solution vector
 */
void qc_socp2inv_pos(pfloat * x, inv_pos_vars * vars, const inv_pos_dims * dims);

#ifdef __cplusplus
}
#endif
  
#endif // __INV_POS_H__
