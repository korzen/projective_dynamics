#include <iostream>
#include <string>
#include <algorithm>
#include <cstring>

extern "C" {
        #include <pk/pk_io.h>
}

#include "pd_mesh.h"

bool arg_flag(char **beg, char **end, const std::string &f){
        return std::find(beg, end, f) != end;
}
// We compile this as C++, so no worries :P
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

int main(int argc, char **argv){
        if (argc < 3 || arg_flag(argv, argv + argc, "-h")){
                printf("Usage: ./pd_mesh_binary <input.json> <output_name> [options]\n"
                                "\t<input.json>     Name of input JSON mesh file\n"
                                "\t<output_name>    Name of output file, binary file will be <output_name>.bmesh\n"
                                "\t-h               Print this help\n");
                return 0;
        }
        char *str = pk_io_read_file(argv[1]);
        assert(str);
        PdMeshSurface *mesh = pd_mesh_surface_mk_from_json(str);
        // Add an attachment if there are none
        if (mesh->n_attachments == 0){
                printf("No attachments found, adding one to attach first index at {0, 0, 0}\n");
                mesh->n_attachments = 1;
                mesh->attachments = (struct PdConstraintAttachment *)malloc(mesh->n_attachments*sizeof *mesh->attachments);
                mesh->attachments[0].i = 0;
                for (size_t i = 0; i < 3; ++i){
                        mesh->attachments[0].position[i] = 0;
                }
        }
        free(str);
        assert(mesh);
        pd_mesh_print_info(mesh);
        std::string out_name = std::string(argv[2]) + ".bmesh";
        if (!pd_mesh_surface_write_binary(mesh, out_name.c_str())){
                printf("Failed to write binary mesh file\n");
                return 1;
        }
        printf("Wrote binary mesh to %s\n", out_name.c_str());
        return 0;
}

