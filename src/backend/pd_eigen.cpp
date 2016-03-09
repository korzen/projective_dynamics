#include <cstdlib>
#include <cstring>

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../pd_time.h"

#include "../pd_solver.h"


struct PdSolver {
        Eigen::VectorXf positions;
        Eigen::VectorXf positions_last;
        Eigen::SparseMatrix<float> mass_mat;
        Eigen::SparseMatrix<float> l_mat;
        Eigen::SparseMatrix<float> j_mat;
        Eigen::SimplicialLLT<Eigen::SparseMatrix<float>, Eigen::Upper> llt;

        struct PdConstraintAttachment *attachments;
        uint32_t n_attachments;

        struct PdConstraintSpring *springs;
        uint32_t n_springs;

        float t2;

        /* cumulative moving average for local and global time */
        double   local_cma;
        double   global_cma;
        uint64_t n_iters;
};



/* compute J matrix */



/*
one iteration algorithm:
1. compute D vector by evaluating constraints (local step; parallel)
2. solve system by plugging int the D vector (global step)
*/


struct PdSolver *
pd_solver_alloc(float const                         *positions,
                uint32_t const                       n_positions, /* TODO: confusing as length of positions == 3*n_positions */
                struct PdConstraintAttachment const *attachments,
                uint32_t const                       n_attachments,
                struct PdConstraintSpring const     *springs,
                uint32_t const                       n_springs,
                float const                          timestep) /* no adaptive timestepping */
{
        struct PdSolver *solver = new struct PdSolver;

        solver->positions = Eigen::VectorXf::Map(positions, 3*n_positions);
        solver->positions_last = solver->positions;


        /* initialize mass matrix */
        float const total_mass = 1.0f;

        /* TODO: not sure if bug in Eigen */
        /*solver->mass_mat.diagonal().setConstant(total_mass/n_positions);*/
        std::vector<Eigen::Triplet<float>> triplets;
        for (uint32_t i = 0; i < n_positions; ++i)
                for (int j = 0; j < 3; ++j)
                        triplets.emplace_back(3*i + j, 3*i + j, total_mass/n_positions);
        solver->mass_mat.resize(3*n_positions, 3*n_positions);
        solver->mass_mat.setFromTriplets(std::begin(triplets), std::end(triplets));


        float const stiffness_attachment = 32.0f;
        float const stiffness_spring     = 16.0f;


        /* build L matrix; kAA^T (and Kronecker product to apply it to 3 vector */
        /* A is spring endpoints incidence matrix as weighted Laplacian matrix */
        /* https://en.wikipedia.org/wiki/Laplacian_matrix */
        triplets.clear();

        /* TODO: not really sure why it is diagonal matrix for attachments since it is stiffness*identity */
        for (uint32_t i = 0; i < n_attachments; ++i) {
                struct PdConstraintAttachment const c = attachments[i];
                for (int j = 0; j < 3; ++j)
                        triplets.emplace_back(3*c.i + j, 3*c.i + j, stiffness_attachment);
        }
        size_t const attachments_size = n_attachments*sizeof *solver->attachments;
        solver->attachments = (struct PdConstraintAttachment *)malloc(attachments_size);
        memcpy(solver->attachments, attachments, attachments_size);
        solver->n_attachments = n_attachments;

        /* have to accumulate, since some entries were filled previously */
        for (uint32_t i = 0; i < n_springs; ++i) {
                struct PdConstraintSpring const c = springs[i];
                for (int j = 0; j < 3; ++j) {
                        /* degree part of L matrix (D matrix) */
                        triplets.emplace_back(3*c.i[0] + j, 3*c.i[0] + j, stiffness_spring);
                        triplets.emplace_back(3*c.i[1] + j, 3*c.i[1] + j, stiffness_spring);

                        /* -incidence part of L matrix (A matrix) */
                        /* since it is simple digraph these are in fact stiffness*-A
                           where A is incidence matrix */
                        triplets.emplace_back(3*c.i[0] + j, 3*c.i[1] + j, -stiffness_spring);
                        triplets.emplace_back(3*c.i[1] + j, 3*c.i[0] + j, -stiffness_spring);
                }
        }
        size_t const springs_size = n_springs*sizeof *solver->springs;
        solver->springs = (struct PdConstraintSpring *)malloc(springs_size);
        memcpy(solver->springs, springs, springs_size);
        solver->n_springs = n_springs;

        solver->l_mat.resize(3*n_positions, 3*n_positions);
        solver->l_mat.setFromTriplets(std::begin(triplets), std::end(triplets));


        /* build J matrix; TODO: explain what the matrix captures */
        /* from paper it is equation 12 */
        /* "ORDER" MUST BE PRESERVED */
        triplets.clear();

        uint32_t offset = 0;
        for (uint32_t i = 0; i < n_attachments; ++i, ++offset) {
                struct PdConstraintAttachment const c = attachments[i];
                for (int j = 0; j < 3; ++j)
                        triplets.emplace_back(3*c.i + j, 3*offset + j, stiffness_attachment);
        }

        for (uint32_t i = 0; i < n_springs; ++i, ++offset) {
                struct PdConstraintSpring const c = springs[i];
                for (int e = 0; e < 2; ++e) {
                        float const k = e ? stiffness_spring : -stiffness_spring;
                        for (int j = 0; j < 3; ++j)
                                triplets.emplace_back(3*c.i[e] + j, 3*offset + j, k);
                }
        }
        uint32_t const n_constraints = n_attachments + n_springs;
        solver->j_mat.resize(3*n_positions, 3*n_constraints);
        solver->j_mat.setFromTriplets(std::begin(triplets), std::end(triplets));

        solver->t2 = timestep*timestep;


        /* prefactor the constant term (mesh topology/masses are fixed) */
        solver->llt.compute(solver->mass_mat + solver->t2*solver->l_mat);


        solver->local_cma  = 0.0;
        solver->global_cma = 0.0;
        solver->n_iters    = 0;


        return solver;
}


void
pd_solver_free(struct PdSolver *solver)
{
        free(solver->attachments);
        free(solver->springs);
        delete solver;
}


void
pd_solver_advance(struct PdSolver *solver)
{
        /* LOCAL STEP (we account everything except the global solve */
        struct timespec local_start;
        clock_gettime(CLOCK_MONOTONIC, &local_start);

        float const gravity = -9.81f;

        /* compute external force */
        Eigen::VectorXf ext_accel = Eigen::VectorXf::Zero(solver->positions.size());
        for (uint32_t i = 0; i < solver->positions.size()/3; ++i)
                ext_accel[3*i + 2] = gravity;
        Eigen::VectorXf const ext_force = solver->mass_mat*ext_accel;


        /* compute d vector for constraints */
        /* MUST GO IN THIS ORDER */
        uint32_t const n_constraints = solver->n_attachments + solver->n_springs;
        Eigen::VectorXf d(3*n_constraints);
        uint32_t offset = 0;
        for (uint32_t i = 0; i < solver->n_attachments; ++i, ++offset) {
                /* TODO: constness */
                struct PdConstraintAttachment c = solver->attachments[i];
                d.block<3, 1>(3*offset, 0) = Eigen::Map<Eigen::Vector3f>(c.position);
        }

        for (uint32_t i = 0; i < solver->n_springs; ++i, ++offset) {
                struct PdConstraintSpring const c = solver->springs[i];
                /* TODO: why swapping indices tamed the crazy net? */
                Eigen::Vector3f const v = solver->positions.block<3, 1>(3*c.i[1], 0) - solver->positions.block<3, 1>(3*c.i[0], 0);
                d.block<3, 1>(3*offset, 0) = v.normalized()*c.rest_length;
        }


        /* solve the global system */
        /* TODO: we can shuffle things to possibly save time, for example add
                 inertia and external force*timestep^2 and the multiply by mass
                 matrix */
        Eigen::VectorXf const y = 2.0f*solver->positions - solver->positions_last;
        Eigen::VectorXf const b = solver->mass_mat*y + solver->t2*(solver->j_mat*d + ext_force);

        /* TODO: this should not copy anything */
        solver->positions_last = solver->positions;

        struct timespec local_end;
        clock_gettime(CLOCK_MONOTONIC, &local_end);


        /* GLOBAL STEP */
        struct timespec global_start;
        clock_gettime(CLOCK_MONOTONIC, &global_start);
        solver->positions = solver->llt.solve(b);

        struct timespec global_end;
        clock_gettime(CLOCK_MONOTONIC, &global_end);

        solver->local_cma = (pd_time_diff_ms(&local_start, &local_end) + solver->n_iters*solver->local_cma)/(solver->n_iters + 1);
        solver->global_cma = (pd_time_diff_ms(&global_start, &global_end) + solver->n_iters*solver->global_cma)/(solver->n_iters + 1);

        if (!(solver->n_iters % 1000)) {
                printf("Local CMA: %f ms\n", solver->local_cma);
                printf("Global CMA: %f ms\n\n", solver->global_cma);
        }

        ++solver->n_iters;
}


float const *
pd_solver_map_positions(struct PdSolver const *solver)
{
        return solver->positions.data();
}
