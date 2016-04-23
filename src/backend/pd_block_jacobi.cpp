#include <cinttypes>
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
        Eigen::SparseMatrix<float> a_mat;
        Eigen::SparseMatrix<float> mass_mat;
        Eigen::SparseMatrix<float> l_mat;
        Eigen::SparseMatrix<float> j_mat;

        struct PdConstraintAttachment *attachments;
        uint32_t n_attachments;

        struct PdConstraintSpring *springs;
        uint32_t n_springs;

        float t2;
        float ext_force[3];

        double   global_time;
        double   local_time;

        /* cumulative moving average for local and global time */
        double   local_cma;
        double   global_cma;
        uint64_t n_iters;

        uint32_t m, n;
        uint32_t block_m, block_n;
};
                

struct PdSolver *
pd_solver_alloc(float const                         *positions,
                uint32_t const                       n_positions,
                struct PdConstraintAttachment const *attachments,
                uint32_t const                       n_attachments,
                struct PdConstraintSpring const     *springs,
                uint32_t const                       n_springs,
                float const                          timestep)
{
        struct PdSolver *solver = new struct PdSolver;

        /* determine number of blocks and their size */
        solver->m = 3;
        solver->n = 3;
        solver->block_m = n_positions;
        solver->block_n = n_positions;
        
        solver->positions = Eigen::VectorXf::Map(positions, 3*n_positions);
        solver->positions_last = solver->positions;

        solver->ext_force[0] = 0.0f;
        solver->ext_force[1] = 0.0f;
        solver->ext_force[2] = -9.8f;

        std::vector<Eigen::Triplet<float>> triplets;

        /* initialize mass matrix */
        float const total_mass = 1.0f;

        /* TODO: not sure if bug in Eigen */
        /*solver->mass_mat.diagonal().setConstant(total_mass/n_positions);*/
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
        solver->a_mat = solver->mass_mat + solver->t2*solver->l_mat;


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
pd_solver_set_ext_force(struct PdSolver *solver, float const *force)
{
        memcpy(solver->ext_force, force, 3*sizeof *force);
}


static void
solve(struct PdSolver *solver, Eigen::VectorXf const &b)
{
        float const epsilon = 1e-1;

        float error = 1.0f;
        uint32_t iters = 0;
        while (error > epsilon) {
                for (uint32_t i = 0; i < solver->m; ++i) {
                        Eigen::VectorXf y = b.block(i*solver->block_n, 0, solver->block_n, 1);

                        for (uint32_t j = 0; j < solver->n; ++j) {
                                if (i != j); else continue;
                                y -= solver->a_mat.block(i*solver->block_m, j*solver->block_n, solver->block_m, solver->block_n)*solver->positions_last.block(j*solver->block_n, 0, solver->block_n, 1);
                        }

                        Eigen::SimplicialLLT<Eigen::SparseMatrix<float>, Eigen::Upper> llt;
                        llt.compute(solver->a_mat.block(i*solver->block_m, i*solver->block_n, solver->block_m, solver->block_n));
                        solver->positions.block(i*solver->block_n, 0, solver->block_n, 1) = llt.solve(y);

                }
                /* TODO: compute error */
                error -= 0.01f;
                ++iters;
        }
        printf("Iteration count: %"PRIu32"\n", iters);
}


void
pd_solver_advance(struct PdSolver *solver, uint32_t const n_iterations)
{
        /* set external force */
        Eigen::VectorXf ext_accel = Eigen::VectorXf::Zero(solver->positions.size());
#pragma omp parallel for
        for (uint32_t i = 0; i < solver->positions.size()/3; ++i)
                for (int j = 0; j < 3; ++j)
                        ext_accel[3*i + j] = solver->ext_force[j];
        Eigen::VectorXf const ext_force = solver->mass_mat*ext_accel;
        Eigen::VectorXf const mass_y = solver->mass_mat*(2.0f*solver->positions - solver->positions_last);

        solver->positions_last = solver->positions;

        for (uint32_t iter = 0; iter < n_iterations; ++iter) {
                /* LOCAL STEP */
                uint32_t const n_constraints = solver->n_attachments + solver->n_springs;
                Eigen::VectorXf d(3*n_constraints);
                for (uint32_t i = 0; i < solver->n_attachments; ++i) {
                        struct PdConstraintAttachment c = solver->attachments[i];
                        d.block<3, 1>(3*i, 0) = Eigen::Map<Eigen::Vector3f>(c.position);
                }

                uint32_t const offset = solver->n_attachments;
                for (uint32_t i = 0; i < solver->n_springs; ++i) {
                        struct PdConstraintSpring c = solver->springs[i];
                        Eigen::Vector3f const v = solver->positions.block<3, 1>(3*c.i[1], 0) - solver->positions.block<3, 1>(3*c.i[0], 0);
                        d.block<3, 1>(3*(offset + i), 0) = v.normalized()*c.rest_length;
                }

                Eigen::VectorXf const b = mass_y + solver->t2*(solver->j_mat*d + ext_force);


                /* GLOBAL STEP */
                /* TODO: pass references to data instead whole solver */
                solve(solver, b);
        }
}


float const *
pd_solver_map_positions(struct PdSolver const *solver)
{
        return solver->positions.data();
}


double
pd_solver_global_cma(struct PdSolver const *solver)
{
        return solver->global_cma;
}


double
pd_solver_local_cma(struct PdSolver const *solver)
{
        return solver->local_cma;
}


char const *
pd_solver_name(struct PdSolver const *solver)
{
        return "Block Jacobi";
}


double
pd_solver_global_time(struct PdSolver const *solver)
{
        return solver->global_time;
}


double
pd_solver_local_time(struct PdSolver const *solver)
{
        return solver->local_time;
}


void
pd_solver_draw_ui(struct PdSolver *solver)
{
}
