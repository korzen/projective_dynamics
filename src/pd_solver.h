#ifndef PD_SOLVER_H
#define PD_SOLVER_H


#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif


struct PdSolver;


struct PdSolver *
pd_solver_alloc(float const *positions,
                uint32_t const n_positions,
                uint8_t const *indices,
                uint32_t const n_indices);


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
