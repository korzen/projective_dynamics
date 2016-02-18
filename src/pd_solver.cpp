#include <cassert>
#include <cstdlib>
#include <cstring>
#include <ctgmath>

#include "pd_solver.h"


static float const gravity = 9.81f;


struct PdSolver {
        float *forces;
        float *positions;
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
        solver->forces = (float *)malloc(3 * n_positions * sizeof *solver->forces);


        /* copy positions */
        solver->positions = (float *)malloc(3 * n_positions * sizeof *positions);
        memcpy(solver->positions, positions, 3 * n_positions * sizeof *positions);
        solver->n_points = n_positions;


        /* determine springs from the mesh edges */
        uint8_t *springs_map = (uint8_t *)calloc(n_indices, sizeof *springs_map);
        for (uint32_t i = 0; i < n_indices; ++i) {
                /* process triangle */
                for (int e = 0; e < 3; ++e) {
                        assert(indices[3 * i + e] < 256);
                        uint8_t spring[2] = {
                                indices[3 * i + e],
                                indices[3 * i  + (e + 1) % 3],
                        };

                        if (spring[0] > spring[1])
                                springs_map[spring[1]] = spring[0];
                        else
                                springs_map[spring[0]] = spring[1];
                }
        }

        /* TODO vector */
        uint32_t n_springs = 0;
        for (uint32_t i = 0; i < n_indices; ++i) {
                if (!springs_map[i])
                        continue;
                ++n_springs;
        }

        solver->springs = (uint8_t *)malloc(2 * n_springs * sizeof *solver->springs);
        solver->n_springs = 0;
        for (uint32_t i = 0; i < n_indices; ++i) {
                if (!springs_map[i])
                        continue;

                solver->springs[2 * solver->n_springs] = i;
                solver->springs[2 * solver->n_springs + 1] = springs_map[i];
                ++solver->n_springs;
        }
        free(springs_map);


        /* compute length of springs */
        solver->lengths = (float *)malloc(solver->n_springs * sizeof *solver->lengths);
        for (uint32_t i = 0; i < solver->n_springs; ++i) {
                float const *p0 = solver->positions + solver->springs[2 * i];
                float const *p1 = solver->positions + solver->springs[2 * i + 1];

                solver->lengths[i] = p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2];
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


void
pd_solver_advance(struct PdSolver *solver)
{
        /* set external point forces */
        for (uint32_t i = 0; i < solver->n_points; ++i) {
                solver->forces[3*i] = 0.0f;
                solver->forces[3*i + 1] = 0.0f;
                solver->forces[3*i + 2] = -gravity;
        }

        /* set spring forces */
        for (uint32_t i = 0; i < solver->n_springs; ++i) {
        }

        /* move points */
}


float const *
pd_solver_map_positions(struct PdSolver const *solver)
{
        return solver->positions;
}
