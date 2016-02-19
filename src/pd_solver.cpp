#include <cassert>
#include <cstdlib>
#include <cstring>
#include <ctgmath>

#include "pd_solver.h"


static float const gravity = 9.81f;


struct PdSolver {
        float *forces;
        float *positions;
        float *positions_last;
        uint32_t n_points;

        uint8_t *springs;
        float *lengths;
        uint32_t n_springs;
};


struct PdSolver *
pd_solver_alloc(float const *positions,
                uint32_t const n_positions,
                uint8_t const *indices,
                uint32_t const n_indices)

{
        struct PdSolver *solver = (struct PdSolver *)malloc(sizeof *solver);
        solver->forces = (float *)malloc(3*n_positions*sizeof *solver->forces);


        /* copy positions */
        solver->positions = (float *)malloc(3*n_positions*sizeof *positions);
        memcpy(solver->positions, positions, 3*n_positions*sizeof *positions);
        solver->positions_last = (float *)malloc(3*n_positions*sizeof *positions);
        memcpy(solver->positions_last, positions, 3*n_positions*sizeof *positions);
        solver->n_points = n_positions;


        /* determine springs from the mesh edges */
        solver->springs = (uint8_t *)malloc(2*n_indices*sizeof *solver->springs);
        solver->n_springs = 0;

        for (uint32_t i = 0; i < n_indices/3; ++i) {
                /* process triangle */
                for (int e = 0; e < 3; ++e) {
                        assert(indices[3*i + e] < 256);
                        uint8_t spring[2] = {
                                indices[3*i + e],
                                indices[3*i  + (e + 1)%3],
                        };

                        if (spring[0] > spring[1]) {
                                uint8_t const tmp = spring[0];
                                spring[0] = spring[1];
                                spring[1] = tmp;
                        }

                        /* skip the spring if the other exists */
                        uint32_t i;
                        for (i = 0; i < solver->n_springs && (solver->springs[2*i] != spring[0] || solver->springs[2*i + 1] != spring[1]); ++i);
                        if (i != solver->n_springs)
                                continue;

                        solver->springs[2*solver->n_springs] = spring[0];
                        solver->springs[2*solver->n_springs + 1] = spring[1];
                        ++solver->n_springs;
                }
        }


        /* compute length of springs */
        solver->lengths = (float *)malloc(solver->n_springs * sizeof *solver->lengths);
        for (uint32_t i = 0; i < solver->n_springs; ++i) {
                float const *p0 = solver->positions + 3*solver->springs[2*i];
                float const *p1 = solver->positions + 3*solver->springs[2*i + 1];
                solver->lengths[i] = sqrt(pow(p1[0] - p0[0], 2) + pow(p1[1] - p0[1], 2) + pow(p1[2] - p0[2], 2));

        }


        return solver;
}


void
pd_solver_free(struct PdSolver *solver)
{
        free(solver->forces);
        free(solver->positions);
        free(solver->springs);
        free(solver->lengths);
        free(solver);
}


/* first vertex is fixed */
void
pd_solver_advance(struct PdSolver *solver, float const timestep)
{
        /* set external point forces */
        for (uint32_t i = 0; i < solver->n_points; ++i) {
                solver->forces[3*i] = 0.0f;
                solver->forces[3*i + 1] = 0.0f;
                solver->forces[3*i + 2] = (i == 0) ? 0.0f : -gravity;
        }

        /* move points with Verlet integration scheme */
        for (uint32_t i = 0; i < 3*solver->n_points; ++i) {
                float const tmp = solver->positions[i];
                /* 2x - x' + a*t^2; x' is last position; assume mass of 1 */
                solver->positions[i] += tmp - solver->positions_last[i] + solver->forces[i] * timestep * timestep;
                solver->positions_last[i] = tmp;
        }

        /* satisfy spring constraints; fixed length assumed; one iteration */
        for (uint32_t i = 0; i < solver->n_springs; ++i) {
                uint8_t const *spring = solver->springs + 2*i;
                float const v[3] = {
                        solver->positions[3*spring[0]] - solver->positions[3*spring[1]],
                        solver->positions[3*spring[0] + 1] - solver->positions[3*spring[1] + 1],
                        solver->positions[3*spring[0] + 2] - solver->positions[3*spring[1] + 2],
                };
                float const v_length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
                float const a = 0.5f*(v_length - solver->lengths[i])/v_length;

                for (int j = 0; j < 3; ++j) {
                        solver->positions[3*spring[0] + j] -= (spring[0] == 0) ? 0.0f : a*v[j];
                        solver->positions[3*spring[1] + j] += (spring[1] == 0) ? 0.0f : a*v[j];
                }

        }
}


float const *
pd_solver_map_positions(struct PdSolver const *solver)
{
        return solver->positions;
}
