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
pd_solver_alloc(float const                         *positions,
                uint32_t const                       n_positions,
                struct PdConstraintAttachment const *attachment_constraints,
                uint32_t const                       n_attachment_constraints,
                struct PdConstraintSpring const     *spring_constraints,
                uint32_t const                       n_spring_constraints,
                float const                          timestep);


void
pd_solver_free(struct PdSolver *solver);


void
pd_solver_set_ext_force(struct PdSolver *solver, float const *force);


void
pd_solver_advance(struct PdSolver *solver);


float const *
pd_solver_map_positions(struct PdSolver const *solver);


/* TODO: should we also provide last frame time? */
double
pd_solver_global_cma(struct PdSolver const *solver);


double
pd_solver_local_cma(struct PdSolver const *solver);


char const *
pd_solver_name(struct PdSolver const *solver);


double
pd_solver_global_time(struct PdSolver const *solver);


double
pd_solver_local_time(struct PdSolver const *solver);

void
pd_solver_draw_ui(struct PdSolver *solver);

#ifdef __cplusplus
}
#endif


#endif
