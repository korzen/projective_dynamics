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


        struct PdConstraintAttachment *attachments;
        uint32_t n_attachments;

        struct PdConstraintSpring *springs;
        uint32_t n_springs;
};


struct PdSolver *
pd_solver_alloc(float const                   *positions,
                uint32_t const                 n_positions,
                struct PdConstraintAttachment *attachment_constraints,
                uint32_t const                 n_attachment_constraints,
                struct PdConstraintSpring     *spring_constraints,
                uint32_t const                 n_spring_constraints)
{
        struct PdSolver *solver = (struct PdSolver *)malloc(sizeof *solver);
        solver->forces = (float *)malloc(3*n_positions*sizeof *solver->forces);

        solver->positions = (float *)malloc(3*n_positions*sizeof *positions);
        memcpy(solver->positions, positions, 3*n_positions*sizeof *positions);
        solver->positions_last = (float *)malloc(3*n_positions*sizeof *positions);
        memcpy(solver->positions_last, positions, 3*n_positions*sizeof *positions);
        solver->n_points = n_positions;

        size_t const attachments_size = n_attachment_constraints*sizeof *solver->attachments;
        solver->attachments = (struct PdConstraintAttachment *)malloc(attachments_size);
        memcpy(solver->attachments, attachment_constraints, attachments_size);
        solver->n_attachments = n_attachment_constraints;

        size_t const springs_size = n_spring_constraints*sizeof *solver->springs;
        solver->springs = (struct PdConstraintSpring *)malloc(springs_size);
        memcpy(solver->springs, spring_constraints, springs_size);
        solver->n_springs = n_spring_constraints;

        return solver;
}


void
pd_solver_free(struct PdSolver *solver)
{
        free(solver->forces);
        free(solver->positions);
        free(solver->attachments);
        free(solver->springs);
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
                solver->positions[i] += tmp - solver->positions_last[i] + solver->forces[i]*timestep*timestep;
                solver->positions_last[i] = tmp;
        }

        /* satisfy spring constraints; fixed length assumed; one iteration */
        for (uint32_t i = 0; i < solver->n_springs; ++i) {
                /* make copy since it is tiny and may get optimized away */
                struct PdConstraintSpring s = solver->springs[i];
                float const v[3] = {
                        solver->positions[3*s.i[0]] - solver->positions[3*s.i[1]],
                        solver->positions[3*s.i[0] + 1] - solver->positions[3*s.i[1] + 1],
                        solver->positions[3*s.i[0] + 2] - solver->positions[3*s.i[1] + 2],
                };
                float const v_length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
                float const a = 0.5f*(v_length - s.rest_length)/v_length;
//                printf("%u, %u %u %f\n", i, s.i[0], s.i[1], s.rest_length);

                for (int j = 0; j < 3; ++j) {
                        solver->positions[3*s.i[0] + j] -= (s.i[0] == 0) ? 0.0f : a*v[j];
                        solver->positions[3*s.i[1] + j] += (s.i[1] == 0) ? 0.0f : a*v[j];
                }

        }
}


float const *
pd_solver_map_positions(struct PdSolver const *solver)
{
        return solver->positions;
}
