#define USE_CUSTOM_KERNELS 1
#define USE_INITIAL_GUESS_Y 0

/* picks first with 1 set */
#define PRECONDITIONER_ILUT        0
#define PRECONDITIONER_AMG         0
#define PRECONDITIONER_JACOBI      1
#define PRECONDITIONER_ROW_SCALING 0
#define PRECONDITIONER_ICHOL0      0
#define PRECONDITIONER_CHOW_PATEL  0
#define PRECONDITIONER_ILU0        0
#define PRECONDITIONER_BLOCK_ILU   0
#define USE_PRECONDITIONER (PRECONDITIONER_ILUT || \
                            PRECONDITIONER_AMG || \
                            PRECONDITIONER_JACOBI || \
                            PRECONDITIONER_ROW_SCALING || \
                            PRECONDITIONER_ICHOL0 || \
                            PRECONDITIONER_CHOW_PATEL || \
                            PRECONDITIONER_ILU0 || \
                            PRECONDITIONER_BLOCK_ILU)

#include <cassert>
#include <vector>
#include <algorithm>
#include <map>

#if USE_CUSPARSE
#include <cusolverSp.h>
#if USE_CUSPARSE_LOW_LEVEL
#include <cusolverSp_LOWLEVEL_PREVIEW.h>
#endif
#endif

#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/amg.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/ichol.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>
#include <imgui/imgui.h>
#include "../pd_time.h"
#include "../pd_solver.h"

#if USE_CUSTOM_KERNELS
#define CU_DEVICE_EXPORT __device__ __host__
#else
#define CU_DEVICE_EXPORT
#endif

// ViennaCL doesn't have small host vectors so we need a lightweight
// vec3f type to do some constraint computation
struct vec3f {
        float x, y, z;

        vec3f(){ vec3f(0.0f, 0.0f, 0.0f); }
        CU_DEVICE_EXPORT vec3f(float x, float y, float z) : x(x), y(y), z(z){}
        vec3f(float x[3]) : x(x[0]), y(x[1]), z(x[2]){}
        float length() const {
                return std::sqrt(x * x + y * y + z * z);
        }
        vec3f normalized() const {
                float len = 1.0f / length();
                return vec3f(x * len, y * len, z * len);
        }
        const float& operator[](size_t i) const {
                switch (i){
                        case 0: return x;
                        case 1: return y;
                        case 2: return z;
                        default: assert(false);
                }
                // Ssh compiler..
                return x;
        }
};
vec3f operator-(const vec3f &a, const vec3f &b){
        return vec3f(a.x - b.x, a.y - b.y, a.z - b.z);
}
vec3f operator*(const vec3f &a, const float c){
        return vec3f(a.x * c, a.y * c, a.z * c);
}

struct PdSolver {
        viennacl::vector<float> positions;
        viennacl::vector<float> positions_last;
        viennacl::compressed_matrix<float> mass_mat;
        viennacl::compressed_matrix<float> l_mat;
        viennacl::compressed_matrix<float> j_mat;
        viennacl::compressed_matrix<float> a_mat;
        
#if PRECONDITIONER_ILUT
        viennacl::linalg::ilut_precond<viennacl::compressed_matrix<float>> *precond_mat;
#elif PRECONDITIONER_AMG
        viennacl::linalg::amg_precond<viennacl::compressed_matrix<float>> *precond_mat;
#elif PRECONDITIONER_JACOBI
        viennacl::linalg::jacobi_precond<viennacl::compressed_matrix<float>> *precond_mat;
#elif PRECONDITIONER_ROW_SCALING
        viennacl::linalg::row_scaling<viennacl::compressed_matrix<float>> *precond_mat;
#elif PRECONDITIONER_ICHOL0
        viennacl::linalg::ichol0_precond<viennacl::compressed_matrix<float>> *precond_mat;
#elif PRECONDITIONER_CHOW_PATEL
        viennacl::linalg::chow_patel_ilu_precond<viennacl::compressed_matrix<float>> *precond_mat;
#elif PRECONDITIONER_ILU0
        viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<float>> *precond_mat;
#elif PRECONDITIONER_BLOCK_ILU
        viennacl::linalg::block_ilu_precond<viennacl::compressed_matrix<float>, viennacl::linalg::ilu0_tag> *precond_mat;
#endif

        std::vector<PdConstraintAttachment> attachments;
        viennacl::vector<float> attachment_pos;
        std::vector<PdConstraintSpring> springs;
        viennacl::vector<int> spring_indices;
        viennacl::vector<float> spring_rest_lengths;

#if USE_CUSPARSE
        cusolverSpHandle_t cusolver_handle;
        cusparseMatDescr_t cusolver_mat_descr;
#if USE_CUSPARSE_LOW_LEVEL
        csrcholInfo_t cusolver_chol_info;
        void *cu_workspace;
        size_t a_mat_size1;
#endif
#else
        float cg_tolerance;
        int cg_max_iterations;
        float cg_last_error;
        int cg_last_iterations;
#endif

        float t2;
        vec3f ext_force;

        // Note: every dereference of a viennacl vector element/iterator
        // triggers a transfer back from the GPU. When mapping we do a single
        // copy into this vector and return a pointer
        std::vector<float> mapped_positions;

        // Most recent global/local times
        double   global_time;
        double   local_time;
        /* cumulative moving average for local and global time */
        double   local_cma;
        double   global_cma;
        uint64_t n_iters;

        std::string name;
};

/* all arrays are copied */
struct PdSolver *
pd_solver_alloc(float const                         *positions,
                uint32_t const                       n_positions,
                struct PdConstraintAttachment const *attachment_constraints,
                uint32_t const                       n_attachment_constraints,
                struct PdConstraintSpring const     *spring_constraints,
                uint32_t const                       n_spring_constraints,
                float const                          timestep)
{
        PdSolver *solver = new PdSolver;
        solver->local_cma = 0.0;
        solver->global_cma = 0.0;
        solver->n_iters = 0;
        solver->t2 = timestep * timestep;
        solver->ext_force = vec3f(0.0f, 0.0f, -9.83f);
#if !USE_CUSPARSE
        solver->cg_tolerance = 1e-5;
        solver->cg_max_iterations = 300;
#endif

        solver->positions = viennacl::vector<float>(3 * n_positions);
        viennacl::copy(positions, positions + 3 * n_positions, solver->positions.begin());
        solver->positions_last = solver->positions;
        solver->mapped_positions.reserve(3 * n_positions);

        // Build the mass matrix, the matrix is a 3n x 3n diagonal matrix
        float const total_mass = 1.0f;
        float const point_mass = total_mass / n_positions;
        // TODO: we could go through set if we want this to be faster
        // Does viennacl have a built in thing for setting diagonal matrices?
        std::vector<std::map<unsigned int, float>> build_sparse_mat(3 * n_positions, std::map<unsigned int, float>());
        for (size_t i = 0; i < build_sparse_mat.size(); ++i){
                build_sparse_mat[i][i] = point_mass;
        }
        solver->mass_mat = viennacl::compressed_matrix<float>(3 * n_positions, 3 * n_positions, 3 * n_positions);
        viennacl::copy(build_sparse_mat, solver->mass_mat);

        float const stiffness_attachment = 32.0f;
        float const stiffness_spring     = 16.0f;

        solver->attachments.reserve(n_attachment_constraints);
        std::copy(attachment_constraints, attachment_constraints + n_attachment_constraints,
                std::back_inserter(solver->attachments));
        solver->attachment_pos = viennacl::vector<float>(3 * n_attachment_constraints);
        {
                std::vector<float> vf_build;
                vf_build.reserve(3 * n_attachment_constraints);
                for (const auto &c : solver->attachments){
                        vf_build.push_back(c.position[0]);
                        vf_build.push_back(c.position[1]);
                        vf_build.push_back(c.position[2]);
                }
                // TODO: We crash here when there are 0 attachments
                // If we add an attachment to the ant then the program just exits? Why?
                // What's different from chihuahua besides size?
                viennacl::copy(vf_build.begin(), vf_build.end(), solver->attachment_pos.begin());
        }

        solver->springs.reserve(n_spring_constraints);
        std::copy(spring_constraints, spring_constraints + n_spring_constraints,
                std::back_inserter(solver->springs));
        solver->spring_indices = viennacl::vector<int>(2 * n_spring_constraints);
        solver->spring_rest_lengths = viennacl::vector<float>(n_spring_constraints);
        {
                std::vector<int> vi_build;
                vi_build.reserve(2 * n_spring_constraints);
                std::vector<float> vf_build;
                vf_build.reserve(n_spring_constraints);
                for (const auto &c : solver->springs){
                        vi_build.push_back(c.i[0]);
                        vi_build.push_back(c.i[1]);
                        vf_build.push_back(c.rest_length);
                }
                viennacl::copy(vi_build.begin(), vi_build.end(), solver->spring_indices.begin());
                viennacl::copy(vf_build.begin(), vf_build.end(), solver->spring_rest_lengths.begin());
        }

        // This is actually building and solving the fast mass spring system not the projective
        // dynamics system. How different are they?
        /* build L matrix; kAA^T (and Kronecker product to apply it to 3 vector */
        /* A is spring endpoints incidence matrix as weighted Laplacian matrix */
        /* https://en.wikipedia.org/wiki/Laplacian_matrix */
        build_sparse_mat.clear();
        build_sparse_mat.resize(3 * n_positions, std::map<unsigned int, float>());
        for (const auto &c : solver->attachments){
                for (size_t j = 0; j < 3; ++j){
                        build_sparse_mat[3 * c.i + j][3 * c.i + j] = stiffness_attachment;
                }
        }
        for (const auto &c : solver->springs){
                for (size_t j = 0; j < 3; ++j){
                        build_sparse_mat[3 * c.i[0] + j][3 * c.i[0] + j] += stiffness_spring;
                        build_sparse_mat[3 * c.i[1] + j][3 * c.i[1] + j] += stiffness_spring;

                        build_sparse_mat[3 * c.i[0] + j][3 * c.i[1] + j] -= stiffness_spring;
                        build_sparse_mat[3 * c.i[1] + j][3 * c.i[0] + j] -= stiffness_spring;
                }
        }
        // Scale all attachments by time squared
        for (auto &r : build_sparse_mat){
                for (auto &v : r){
                        v.second *= solver->t2;
                }
        }
        solver->l_mat = viennacl::compressed_matrix<float>(3 * n_positions, 3 * n_positions);
        viennacl::copy(build_sparse_mat, solver->l_mat);
        // Add in mass to form the A matrix for the system
        for (size_t i = 0; i < build_sparse_mat.size(); ++i){
                build_sparse_mat[i][i] += point_mass;
        }
        solver->a_mat = viennacl::compressed_matrix<float>(3 * n_positions, 3 * n_positions);
        viennacl::copy(build_sparse_mat, solver->a_mat);

        // Build the J matrix
        // TODO: I'm not sure what pavol's code means here
        build_sparse_mat.clear();
        build_sparse_mat.resize(3 * n_positions, std::map<unsigned int, float>());
        // TODO what is offset for in pavol's code?
        // the offset guarantess that springs go right after attachments (could use solver->n_attachments)
        size_t offset = 0;
        for (const auto &c : solver->attachments){
                for (size_t j = 0; j < 3; ++j){
                    build_sparse_mat[3 * c.i + j][3 * offset + j] = stiffness_attachment;
                }
                ++offset;
        }
        for (const auto &c : solver->springs){
                for (size_t e = 0; e < 2; ++e){
                        float const k = e ? stiffness_spring : -stiffness_spring;
                        for (size_t j = 0; j < 3; ++j){
                            build_sparse_mat[3 * c.i[e] + j][3 * offset + j] += k;
                        }
                }
                ++offset;
        }
        solver->j_mat = viennacl::compressed_matrix<float>(3 * n_positions, 3 * offset);
        viennacl::copy(build_sparse_mat, solver->j_mat);

#if !USE_CUSPARSE

#ifdef VIENNACL_WITH_CUDA
        solver->name = "ViennaCL with CUDA";
#elif defined(VIENNACL_WITH_OPENCL)
        solver->name = "ViennaCL with OpenCL";
#else
        solver->name = "ViennaCL";
#endif

#else

#if !USE_CUSPARSE_LOW_LEVEL
        solver->name = "cuSPARSE High Level API";
#else
        solver->name = "cuSPARSE Low Level API";
#endif

#endif

        // TODO: preconditioner settings
#if PRECONDITIONER_ILUT
        const viennacl::linalg::ilut_tag ilut_conf(100, 1e-10);
        solver->precond_mat = new viennacl::linalg::ilut_precond<viennacl::compressed_matrix<float>>(solver->a_mat, ilut_conf);
        solver->name += " precond ILUT";
#elif PRECONDITIONER_AMG
        viennacl::linalg::amg_tag amg_conf;
        solver->precond_mat = new viennacl::linalg::amg_precond<viennacl::compressed_matrix<float>>(solver->a_mat, amg_conf);
        solver->precond_mat->setup();
        solver->name += " precond AMG";
#elif PRECONDITIONER_JACOBI
        solver->precond_mat = new viennacl::linalg::jacobi_precond<viennacl::compressed_matrix<float>>(solver->a_mat,
                                                                                                       viennacl::linalg::jacobi_tag());
        solver->name += " precond Jacobi";
#elif PRECONDITIONER_ROW_SCALING
        solver->precond_mat = new viennacl::linalg::row_scaling<viennacl::compressed_matrix<float>>(solver->a_mat,
                                                                                                    viennacl::linalg::row_scaling_tag());
        solver->name += " precond Row Scaling";
#elif PRECONDITIONER_ICHOL0
        solver->precond_mat = new viennacl::linalg::ichol0_precond<viennacl::compressed_matrix<float>>(solver->a_mat,
                                                                                                       viennacl::linalg::ichol0_tag());
        solver->name += " precond ICHOL0";
#elif PRECONDITIONER_CHOW_PATEL
        /* first argument is number of nonlinear sweeps; second is number of Jacobi iterations per triangular "solve" */
        viennacl::linalg::chow_patel_tag chow_patel_conf(3, 2);
        solver->precond_mat = new viennacl::linalg::chow_patel_ilu_precond<viennacl::compressed_matrix<float>>(solver->a_mat,
                                                                                                               chow_patel_conf);
        solver->name += " precond Chow Patel";
#elif PRECONDITIONER_ILU0
        viennacl::linalg::ilu0_tag ilu0_conf;
        solver->precond_mat = new viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<float>>(solver->a_mat, ilu0_conf);
        solver->name += " precond ILU0";
#elif PRECONDITIONER_BLOCK_ILU
        viennacl::linalg::ilu0_tag ilu0_conf;
        solver->precond_mat = new viennacl::linalg::block_ilu_precond<viennacl::compressed_matrix<float>, viennacl::linalg::ilu0_tag>(solver->a_mat,
                                                                                                                                      ilu0_conf);
        solver->name += " precond Block ILU";
#endif

        std::cout << "NNZ in A matrix: " << solver->a_mat.nnz() << '\n';
        std::cout << "Sparsity of A matrix: " << (double)solver->a_mat.nnz()/solver->a_mat.size1()/solver->a_mat.size2() << '\n';

#if USE_CUSPARSE
        std::cout << "Setting up for cuSPARSE solver\n";
        // Setup cuSparse solver
        if (cusolverSpCreate(&solver->cusolver_handle) != CUSOLVER_STATUS_SUCCESS){
                std::cout << "cuda error creating solver\n";
        }
        if (cusparseCreateMatDescr(&solver->cusolver_mat_descr) != CUSPARSE_STATUS_SUCCESS){
                std::cout << "cuda error creating mat description\n";
        }
        // TODO: We say general here, but should really say symmetric and then only
        // keep the upper or lower triangle of A (general is the default)
        //cusparseSetMatType(solver->cusolver_mat_descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC);
        cusparseSetMatDiagType(solver->cusolver_mat_descr, CUSPARSE_DIAG_TYPE_NON_UNIT);

#if USE_CUSPARSE_LOW_LEVEL
        std::cout << "Setting up for low level cuSPARSE solver\n";
        if (cusolverSpCreateCsrcholInfo(&solver->cusolver_chol_info) != CUSOLVER_STATUS_SUCCESS){
                std::cout << "Error creating cusolver chol info\n";
        }
        // Perform symoblic analysis on the matrix
        auto status = cusolverSpXcsrcholAnalysis(solver->cusolver_handle, solver->a_mat.size1(),
                solver->a_mat.nnz(), solver->cusolver_mat_descr,
                viennacl::cuda_arg<int>(solver->a_mat.handle1()),
                viennacl::cuda_arg<int>(solver->a_mat.handle2()), solver->cusolver_chol_info);
        if (status != CUSOLVER_STATUS_SUCCESS){
                std::cout << "error performing symbolic analysis\n";
        }
        size_t cu_internal_data_in_bytes = 0;
        size_t cu_workspace_in_bytes = 0;
        // Setup room for the factorization working space
        status = cusolverSpScsrcholBufferInfo(solver->cusolver_handle, solver->a_mat.size1(),
                solver->a_mat.nnz(), solver->cusolver_mat_descr,
                viennacl::cuda_arg<float>(solver->a_mat.handle()),
                viennacl::cuda_arg<int>(solver->a_mat.handle1()),
                viennacl::cuda_arg<int>(solver->a_mat.handle2()),
                solver->cusolver_chol_info, &cu_internal_data_in_bytes, &cu_workspace_in_bytes);
        if (status != CUSOLVER_STATUS_SUCCESS){
                std::cout << "error performing cuda buffer info\n";
        }
        if (cudaMalloc(&solver->cu_workspace, cu_workspace_in_bytes) != cudaSuccess){
                std::cout << "error allocating cuda mem\n";
        }
        // Next thou shalt factorize
        status = cusolverSpScsrcholFactor(solver->cusolver_handle, solver->a_mat.size1(),
                solver->a_mat.nnz(), solver->cusolver_mat_descr,
                viennacl::cuda_arg<float>(solver->a_mat.handle()),
                viennacl::cuda_arg<int>(solver->a_mat.handle1()),
                viennacl::cuda_arg<int>(solver->a_mat.handle2()),
                solver->cusolver_chol_info, solver->cu_workspace);
        if (status != CUSOLVER_STATUS_SUCCESS){
                std::cout << "cusolver error factorizing matrix\n";
        }
        solver->a_mat_size1 = solver->a_mat.size1();
        // Dump the old matrix
        solver->a_mat.resize(1, 1, false);
#endif
#endif
        return solver;
}

void
pd_solver_free(struct PdSolver *solver)
{
#if USE_CUSPARSE
        cusolverSpDestroy(solver->cusolver_handle);
        cusparseDestroyMatDescr(solver->cusolver_mat_descr);
#if USE_CUSPARSE_LOW_LEVEL
        cusolverSpDestroyCsrcholInfo(solver->cusolver_chol_info);
        cudaFree(solver->cu_workspace);
#endif
#endif

#if USE_PRECONDITIONER
        delete solver->precond_mat;
#endif
        delete solver;
}

#if defined(VIENNACL_WITH_CUDA) && USE_CUSTOM_KERNELS
// Kernel to set up the external acceleration vector
__global__ void set_external_acceleration(float *out, const int size, const vec3f ext_force){
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < size){
                out[3 * i] = ext_force.x;
                out[3 * i + 1] = ext_force.y;
                out[3 * i + 2] = ext_force.z;
        }
}
// Kernel to set the position attachment constraints
__global__ void set_attachment_constraints(const float * __restrict__ attachments, float *out, const int size){
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < size){
                out[3 * i] = attachments[3 * i];
                out[3 * i + 1] = attachments[3 * i + 1];
                out[3 * i + 2] = attachments[3 * i + 2];
        }
}
// TODO WILL: These functions aren't provided by CUDA?
__device__ float3 operator-(const float3 a, const float3 b){
        return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
__device__ float3 operator*(const float3 a, const float s){
        return make_float3(a.x * s, a.y * s, a.z * s);
}
__device__ float3 normalize(const float3 a){
        float len = rsqrt(a.x * a.x + a.y * a.y + a.z * a.z);
        return make_float3(a.x * len, a.y * len, a.z * len);
}
// Kernel to set the spring constraints
__global__ void set_spring_constraints(const float * __restrict__ positions, const int * __restrict__ indices,
        const float * __restrict__ lengths, const int offset, float * __restrict__ out, const int size)
{
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < size){
                const int a = indices[2 * i];
                const int b = indices[2 * i + 1];
                const float3 pa = make_float3(positions[3 * a], positions[3 * a + 1], positions[3 * a + 2]);
                const float3 pb = make_float3(positions[3 * b], positions[3 * b + 1], positions[3 * b + 2]);
                const float3 v = normalize(pb - pa) * lengths[i];
                out[3 * (i + offset)] = v.x;
                out[3 * (i + offset) + 1] = v.y;
                out[3 * (i + offset) + 2] = v.z;
        }
}

#endif


void
pd_solver_set_ext_force(struct PdSolver *solver, const float *force)
{
        solver->ext_force = vec3f(force[0], force[1], force[2]);
}


void
pd_solver_advance(struct PdSolver *solver){
        /* LOCAL STEP (we account everything except the global solve */
        struct timespec local_start;
        clock_gettime(CLOCK_MONOTONIC, &local_start);


        // Compute external forces (gravity)
        std::vector<float> vec_build(solver->positions.size(), 0.0f);
        viennacl::vector<float> ext_force(solver->positions.size());
#if !defined(VIENNACL_WITH_CUDA) || !USE_CUSTOM_KERNELS
        // TODO: For OpenCL we need to also write these init kernels
        for (size_t i = 0; i < solver->positions.size() / 3; ++i){
                vec_build[3 * i] = solver->ext_force.x;
                vec_build[3 * i + 1] = solver->ext_force.y;
                vec_build[3 * i + 2] = solver->ext_force.z;
        }
        viennacl::copy(vec_build.begin(), vec_build.end(), ext_force.begin());
#else
        set_external_acceleration<<<(solver->positions.size() / 3) / 32 + 1, 32>>>(viennacl::cuda_arg(ext_force),
                solver->positions.size() / 3, solver->ext_force);
#endif
        ext_force = viennacl::linalg::prod(solver->mass_mat, ext_force);

        // Setup the constraints vector
        size_t const n_constraints = solver->attachments.size() + solver->springs.size();
        viennacl::vector<float> d(3 * n_constraints);
#if !defined(VIENNACL_WITH_CUDA) || !USE_CUSTOM_KERNELS
        vec_build.clear();
        vec_build.resize(3 * n_constraints, 0);
        size_t offset = 0;
        for (const auto &c : solver->attachments){
                for (size_t j = 0; j < 3; ++j){
                        vec_build[3 * offset + j] = c.position[j];
                }
                ++offset;
        }
        // TODO: We really should compute constraints on gpu
        float const *mapped_pos = pd_solver_map_positions(solver);
        for (const auto &c : solver->springs){
                vec3f const pa = vec3f(mapped_pos[3 * c.i[1]], mapped_pos[3 * c.i[1] + 1], mapped_pos[3 * c.i[1] + 2]);
                vec3f const pb = vec3f(mapped_pos[3 * c.i[0]], mapped_pos[3 * c.i[0] + 1], mapped_pos[3 * c.i[0] + 2]);
                vec3f const v = (pa - pb).normalized() * c.rest_length;
                for (size_t j = 0; j < 3; ++j){
                        vec_build[3 * offset + j] = v[j];
                }
                ++offset;
        }
        viennacl::copy(vec_build.begin(), vec_build.end(), d.begin());
#else
        set_attachment_constraints<<<solver->attachments.size() / 32 + 1, 32>>>(viennacl::cuda_arg(solver->attachment_pos),
                    viennacl::cuda_arg(d), solver->attachments.size());

        set_spring_constraints<<<solver->springs.size() / 32 + 1, 32>>>(viennacl::cuda_arg(solver->positions),
                viennacl::cuda_arg(solver->spring_indices), viennacl::cuda_arg(solver->spring_rest_lengths),
                solver->attachments.size(), viennacl::cuda_arg(d), solver->springs.size());
#endif

        viennacl::vector<float> y = 2.0f * solver->positions - solver->positions_last;
        viennacl::vector<float> b = viennacl::linalg::prod(solver->mass_mat, y) + solver->t2
            * (viennacl::linalg::prod(solver->j_mat, d) + ext_force);

        solver->positions_last = solver->positions;

        struct timespec local_end;
        clock_gettime(CLOCK_MONOTONIC, &local_end);

        /* GLOBAL STEP */
        struct timespec global_start;
        clock_gettime(CLOCK_MONOTONIC, &global_start);
        const viennacl::linalg::cg_tag custom_cg(solver->cg_tolerance, solver->cg_max_iterations);
        viennacl::linalg::cg_solver<viennacl::vector<float>> custom_solver(custom_cg);

#if USE_INITIAL_GUESS_Y
        custom_solver.set_initial_guess(y);
#endif

#if !USE_CUSPARSE
        // Solve the system with ViennaCL's CG solver
#if USE_PRECONDITIONER
        solver->positions = custom_solver(solver->a_mat, b, *solver->precond_mat);
#else
        // default cg uses tolerance 1e-8 and at most 300 iterations
        solver->positions = custom_solver(solver->a_mat, b);
#endif
        solver->cg_last_iterations = custom_solver.tag().iters();
        solver->cg_last_error = custom_solver.tag().error();

#else
#if !USE_CUSPARSE_LOW_LEVEL
        // Solve the system with cuSPARSE's higher level API
        int a_mat_singularity = 0;
        auto status = cusolverSpScsrlsvchol(solver->cusolver_handle, solver->a_mat.size1(), solver->a_mat.nnz(),
                solver->cusolver_mat_descr, viennacl::cuda_arg<float>(solver->a_mat.handle()),
                viennacl::cuda_arg<int>(solver->a_mat.handle1()), viennacl::cuda_arg<int>(solver->a_mat.handle2()),
                viennacl::cuda_arg(b), 0.0001, 0,
                viennacl::cuda_arg(solver->positions), &a_mat_singularity);
#else
        // Solve the system with cuSPARSE's low level API
        auto status = cusolverSpScsrcholSolve(solver->cusolver_handle, solver->a_mat_size1,
                viennacl::cuda_arg(b), viennacl::cuda_arg(solver->positions), solver->cusolver_chol_info,
                solver->cu_workspace);
#endif
        // Wait for cuda solve to be done
        cudaDeviceSynchronize();
        if (status != CUSOLVER_STATUS_SUCCESS){
                std::cout << "cuda error solving the global system\n";
        }
#endif
        struct timespec global_end;
        clock_gettime(CLOCK_MONOTONIC, &global_end);

        solver->global_time = pd_time_diff_ms(&global_start, &global_end);
        solver->local_time = pd_time_diff_ms(&local_start, &local_end);

        solver->global_cma = (solver->global_time + solver->n_iters*solver->global_cma)/(solver->n_iters + 1);
        solver->local_cma = (solver->local_time + solver->n_iters*solver->local_cma)/(solver->n_iters + 1);

        if (solver->n_iters && !(solver->n_iters % 500)) {
                printf("Local CMA: %f ms\n", solver->local_cma);
                printf("Global CMA: %f ms\n\n", solver->global_cma);
        }
        ++solver->n_iters;
}

float const *
pd_solver_map_positions(struct PdSolver const *solver){
        PdSolver *sv = const_cast<PdSolver*>(solver);
        sv->mapped_positions.clear();
        sv->mapped_positions.resize(solver->positions.size(), 0);
        viennacl::copy(solver->positions.begin(), solver->positions.end(), sv->mapped_positions.begin());
        return sv->mapped_positions.data();
}

double
pd_solver_global_cma(const struct PdSolver *solver){
        return solver->global_cma;
}

double
pd_solver_local_cma(const struct PdSolver *solver){
        return solver->local_cma;
}

char const *
pd_solver_name(struct PdSolver const *solver)
{
        return solver->name.c_str();
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
#if !USE_CUSPARSE
        if (ImGui::Begin("GPU Solver")){
                ImGui::Text("Solver: %s", pd_solver_name(solver));
                ImGui::InputFloat("CG Tolerance", &solver->cg_tolerance, 0, 0, 8);

                ImGui::InputInt("CG Max Iterations", &solver->cg_max_iterations);
                // Make sure we do at least on CG iteration
                solver->cg_max_iterations = std::max(solver->cg_max_iterations, 1);
                ImGui::Text("Last Solve: %d iterations %f error", solver->cg_last_iterations,
                        solver->cg_last_error);
        }
        ImGui::End();
#endif
}
