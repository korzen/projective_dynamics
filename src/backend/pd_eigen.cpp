#include <cstdlib>
#include <cstring>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../pd_solver.h"


struct PdSolver {
        Eigen::VectorXf positions;
        Eigen::VectorXf positions_last;
        Eigen::SparseMatrix<float> mass_mat;
        Eigen::SparseMatrix<float> l_mat;
        Eigen::SparseMatrix<float> j_mat;

        struct PdConstraintAttachment *attachments;
        uint32_t n_attachments;

        struct PdConstraintSpring *springs;
        uint32_t n_springs;
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
                uint32_t const                       n_springs)
{
        struct PdSolver *solver = new struct PdSolver;

        solver->positions = Eigen::VectorXf::Map(positions, 3*n_positions);
        solver->positions_last = solver->positions;


        /* initialize mass matrix */
        float const total_mass = 1.0f;
        solver->mass_mat.resize(3*n_positions, 3*n_positions);

        /* TODO: not sure if bug in Eigen */
        /*solver->mass_mat.diagonal().setConstant(total_mass/n_positions);*/
        for (uint32_t i = 0; i < n_positions; ++i)
                for (int j = 0; j < 3; ++j)
                solver->mass_mat.coeffRef(3*i + j, 3*i + j) = total_mass/n_positions;


        float const stiffness_attachment = 32.0f;
        float const stiffness_spring     = 16.0f;


        /* build L matrix; kAA^T (and Kronecker product to apply it to 3 vector */
        /* A is spring endpoints incidence matrix as weighted Laplacian matrix */
        /* https://en.wikipedia.org/wiki/Laplacian_matrix */
        solver->l_mat.resize(3*n_positions, 3*n_positions);

        /* TODO: not really sure why it is diagonal matrix for attachments since it is stiffness*identity */
        for (uint32_t i = 0; i < n_attachments; ++i) {
                struct PdConstraintAttachment const c = attachments[i];
                for (int j = 0; j < 3; ++j)
                        solver->l_mat.coeffRef(3*c.i + j, 3*c.i + j) += stiffness_attachment;
        }
        size_t const attachments_size = n_attachments*sizeof *solver->attachments;
        solver->attachments = (struct PdConstraintAttachment *)malloc(attachments_size);
        memcpy(solver->attachments, attachments, attachments_size);
        solver->n_attachments = n_attachments;

        /* have to accumulate, since some entries were filled previously */
        /* TODO: build triplets to avoid potential bugs */
        for (uint32_t i = 0; i < n_springs; ++i) {
                struct PdConstraintSpring const c = springs[i];
                for (int j = 0; j < 3; ++j) {
                        /* degree part of L matrix (D matrix) */
                        solver->l_mat.coeffRef(3*c.i[0] + j, 3*c.i[0] + j) += stiffness_spring;
                        solver->l_mat.coeffRef(3*c.i[1] + j, 3*c.i[1] + j) += stiffness_spring;

                        /* -incidence part of L matrix (A matrix) */
                        /* since it is simple digraph these are in fact stiffness*-A
                           where A is incidence matrix */
                        solver->l_mat.coeffRef(3*c.i[0] + j, 3*c.i[1] + j) -= stiffness_spring;
                        solver->l_mat.coeffRef(3*c.i[1] + j, 3*c.i[0] + j) -= stiffness_spring;
                }
        }
        size_t const springs_size = n_springs*sizeof *solver->springs;
        solver->springs = (struct PdConstraintSpring *)malloc(springs_size);
        memcpy(solver->springs, springs, springs_size);
        solver->n_springs = n_springs;


        /* build J matrix; TODO: explain what the matrix captures */
        /* from paper it is equation 12 */
        /* "ORDER" MUST BE PRESERVED */
        uint32_t const n_constraints = n_attachments + n_springs;
        solver->j_mat.resize(3*n_positions, 3*n_constraints);

        /* TODO: building triples may be faster */
        uint32_t offset = 0;
        for (uint32_t i = 0; i < n_attachments; ++i, ++offset) {
                struct PdConstraintAttachment const c = attachments[i];
                /* TODO: what stiffness for attachment? */
                solver->j_mat.coeffRef(3*c.i, 3*offset) = stiffness_attachment;
                solver->j_mat.coeffRef(3*c.i + 1, 3*offset + 1) = stiffness_attachment;
                solver->j_mat.coeffRef(3*c.i + 2, 3*offset + 2) = stiffness_attachment;
        }

        for (uint32_t i = 0; i < n_springs; ++i, ++offset) {
                struct PdConstraintSpring const c = springs[i];
                for (int e = 0; e < 2; ++e) {
                        float const k = e ? stiffness_spring : -stiffness_spring;
                        solver->j_mat.coeffRef(3*c.i[e], 3*offset) = k;
                        solver->j_mat.coeffRef(3*c.i[e] + 1, 3*offset + 1) = k;
                        solver->j_mat.coeffRef(3*c.i[e] + 2, 3*offset + 2) = k;
                }
        }


        return solver;
}


void
pd_solver_free(struct PdSolver *solver)
{
        delete solver;
}


void
pd_solver_advance(struct PdSolver *solver, float const timestep)
{
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
        Eigen::VectorXf const b = solver->mass_mat*y + timestep*timestep*(solver->j_mat*d + ext_force);

        /* TODO: this should not copy anything */
        solver->positions_last = solver->positions;

        Eigen::SimplicialLLT<Eigen::SparseMatrix<float>, Eigen::Upper> llt_solver;
        llt_solver.compute(solver->mass_mat + timestep*timestep*solver->l_mat);
        solver->positions = llt_solver.solve(b);
}


float const *
pd_solver_map_positions(struct PdSolver const *solver)
{
        return solver->positions.data();
}
