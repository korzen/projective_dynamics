#ifndef PD_SOLVER_H
#define PD_SOLVER_H


#include <stdint.h>

#include "pd_constraint.h"


#ifdef __cplusplus
extern "C" {
#endif


struct PdSolver;


/* all arrays are copied */
struct PdSolver *
pd_solver_alloc(float const                   positions[],
                uint32_t const                n_positions,
                struct PdConstraintAttachment attachment_constraints[],
                uint32_t const                n_attachment_constraints,
                struct PdConstraintSpring     spring_constraints[],
                uint32_t const                n_spring_constraints);


void
pd_solver_free(struct PdSolver *solver);


void
pd_solver_advance(struct PdSolver *solver, float const timestep);


float const *
pd_solver_map_positions(struct PdSolver const *solver);


#ifdef __cplusplus
}
#endif


#endif
