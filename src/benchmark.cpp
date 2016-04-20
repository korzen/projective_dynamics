#include <iostream>
#include <array>
#include <string>
#include <cstring>
#include <algorithm>
#include <ratio>
#include <sstream>

#include <pico_bench/pico_bench.h>

extern "C" {
        #include <pk/pk_io.h>
        #include <pk/pk_linalg.h>
}

#include "pd_mesh.h"
#include "pd_solver.h"

namespace pb = pico_bench;

bool arg_flag(char **beg, char **end, const std::string &f){
        return std::find(beg, end, f) != end;
}
template<typename T>
T get_arg(char **beg, char **end, const std::string &f){
        char **it = std::find(beg, end, f);
        if (it != end && ++it != end){
                std::stringstream ss;
                ss << *it;
                T t;
                ss >> t;
                return t;
        }
        return T();
}
template<typename T, size_t N>
std::array<T, N> get_arg(char **beg, char **end, const std::string &f){
        char **it = std::find(beg, end, f);
        assert(it + 1 != end);
        ++it;
        std::array<T, N> arr;
        for (size_t read = 0; read < N && it != end; ++read, ++it){
                std::stringstream ss;
                ss << *it;
                ss >> arr[read];
        }
        return arr;
}

int main(int argc, char **argv){
        // Check that a valid iterations and time have been passed
        int bench_iters = 0;
        int bench_seconds = 0;
        if (argc >= 3){
                bench_iters = std::atoi(argv[1]);
                bench_seconds = std::atoi(argv[2]);
        }
        if (argc < 3 || bench_iters <= 0 || bench_seconds <= 0){
            std::cout << "Usage: ./pd_benchmark <max iters> <max time in seconds> [options]\n"
                << "\tmax iters and max time in seconds must be >= 0\nOptions:\n"
                << "\t--size <x> <y>       Cloth mesh size\n"
                << "\t--mesh <filename>    Tet mesh file to load\n"
                << "\t-n <number>          Number of iterations of projective dynamics per timestep\n";
            return 1;
        }

        std::array<int, 2> cloth_resolution = {10, 10};
        std::string mesh_filename;
        uint32_t n_iterations = 10;

        if (arg_flag(argv, argv + argc, "--size")){
                cloth_resolution = get_arg<int, 2>(argv, argv + argc, "--size");
                std::cout << "Using cloth of size " << cloth_resolution[0] << "x"
                    << cloth_resolution[1] << "\n";
        } else if (arg_flag(argv, argv + argc, "--mesh")){
                mesh_filename = get_arg<std::string>(argv, argv + argc, "--mesh");
        }
        if (arg_flag(argv, argv + argc, "-n")){
                n_iterations = get_arg<uint32_t>(argv, argv + argc, "-n");
        }

        const float timestep = 1.0f / (60.0f * n_iterations);

        PdMeshSurface *mesh = nullptr;
        if (!mesh_filename.empty()){
                char *str = pk_io_read_file(mesh_filename.c_str());
                assert(str);
                mesh = pd_mesh_surface_mk_from_json(str);
                free(str);
        } else {
                mesh = pd_mesh_surface_mk_grid(cloth_resolution[0], cloth_resolution[1]);
        }

        assert(mesh);
        pd_mesh_print_info(mesh);
        PdSolver *solver = pd_solver_alloc(mesh->positions, mesh->n_positions,
                mesh->attachments, mesh->n_attachments, mesh->springs, mesh->n_springs, timestep);
        assert(solver);

        std::cout << "Benchmarking solver " << pd_solver_name(solver) << "\n";

        using millis = std::chrono::duration<double, std::milli>;
        pb::Benchmarker<millis> bencher(bench_iters, std::chrono::seconds{bench_seconds});
        // We're also custom timing the local and global steps as well
        std::vector<millis> local_steps;
        std::vector<millis> global_steps;
        local_steps.reserve(bench_iters * n_iterations);
        global_steps.reserve(bench_iters * n_iterations);
        auto full_solve_stats = bencher([&](){
                auto start = std::chrono::high_resolution_clock::now();
                for (uint32_t i = 0; i < n_iterations; ++i){
                        pd_solver_advance(solver);

                        local_steps.push_back(millis{pd_solver_local_time(solver)});
                        global_steps.push_back(millis{pd_solver_global_time(solver)});
                }
                auto end = std::chrono::high_resolution_clock::now();
                return std::chrono::duration_cast<millis>(end - start);
        });
        pb::Statistics<millis> local_stats{std::move(local_steps)};
        pb::Statistics<millis> global_stats{std::move(global_steps)};
        local_stats.time_suffix = " ms";
        global_stats.time_suffix = " ms";
        full_solve_stats.time_suffix = " ms";
        pd_solver_free(solver);

        pd_mesh_surface_free(mesh);

        std::cout << "Backend: " << pd_solver_name(solver) << "\n"
                << "Local Solve stats:\n" << local_stats << "\n"
                << "Global Solve stats:\n" << global_stats << "\n"
                << "Full Solve stats:\n" << full_solve_stats << "\n";
        return 0;
}
