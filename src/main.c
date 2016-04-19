#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#include <algorithm>
#include <array>

#include <epoxy/gl.h>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw_gl3.h>
#include <GLFW/glfw3.h>

extern "C" {
        #include <pk/pk_io.h>
        #include <pk/pk_linalg.h>
}

#include "pd_mesh.h"
#include "pd_solver.h"


enum { EBO_TRIANGLES, EBO_LINES, N_EBOS, };

static GLuint  ebos[N_EBOS];
static GLsizei lines_count;
static GLuint  pipeline;
static GLuint  programs[2];
static GLsizei triangles_count;
static GLuint  ubo;
static GLuint  vao;
static GLuint  vbo;

static quat_t cur_quat;
static uint32_t n_positions;
static float *positions_mapped;
static struct MatBlock {
        mat4_t model;
        mat4_t view;
        mat4_t projection;
} *ubo_mapped;

enum { n_solvers = 2, };
static struct PdSolver *solvers[n_solvers];

static uint32_t n_iterations = 10;
static float timestep = 1.0f/(60.0f*10);
static uint32_t resolution_x = 16;
static uint32_t resolution_y = 16;
static char *mesh_filename;

/* view control matrices */
static mat4_t rotation = MAT4;
static mat4_t panning  = MAT4;
static mat4_t zoom     = MAT4;

static bool hover = false;
static vec4_t gravity = VEC4(0.0f, 0.0f, -9.8f, 0.0f);
static double x_pos_prev, y_pos_prev;


static quat_t
mouse_quat(GLFWwindow *window)
{
        double x_win, y_win;
        glfwGetCursorPos(window, &x_win, &y_win);

        int width, height;
        glfwGetWindowSize(window, &width, &height);


        /* convert to NDC */
        float const x = 2.0f*x_win/width - 1.0f;
        float const y = 1.0f - 2.0f*y_win/height;

        float const z2 = 1.0f - x*x - y*y;
        if (z2 < 0.0f) {
                float const length = sqrt(x*x + y*y);
                return QUAT(0.0f, x/length, y/length, 0.0f);
        }

        return QUAT(0.0f, x, y, sqrt(z2));
}


static void
debug(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length,
      GLchar const *message, void const *user_param)
{
        const char *severity_msg = NULL;
        switch (severity){
            case GL_DEBUG_SEVERITY_HIGH_ARB:
                severity_msg = "High severity";
                break;
            case GL_DEBUG_SEVERITY_MEDIUM_ARB:
                severity_msg = "Medium severity";
                break;
            case GL_DEBUG_SEVERITY_LOW_ARB:
                severity_msg = "Low severity";
                break;
            default:
                severity_msg = "Unknown severity";
        }
        const char *src_msg = NULL;
        switch (source){
            case GL_DEBUG_SOURCE_API_ARB:
                src_msg = "API";
                break;
            case GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB:
                src_msg = "Window system";
                break;
            case GL_DEBUG_SOURCE_SHADER_COMPILER_ARB:
                src_msg = "Shader compiler";
                break;
            case GL_DEBUG_SOURCE_THIRD_PARTY_ARB:
                src_msg = "Third party";
                break;
            case GL_DEBUG_SOURCE_APPLICATION_ARB:
                src_msg = "Application";
                break;
            default:
                src_msg = "Other";
        }
        const char *type_msg = NULL;
        switch (type){
            case GL_DEBUG_TYPE_ERROR_ARB:
                type_msg = "Error";
                break;
            case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB:
                type_msg = "Deprecated behavior";
                break;
            case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB:
                type_msg = "Undefined behavior";
                break;
            case GL_DEBUG_TYPE_PORTABILITY_ARB:
                type_msg = "Portability";
                break;
            case GL_DEBUG_TYPE_PERFORMANCE_ARB:
                type_msg = "Performance";
                break;
            default:
                type_msg = "Other";
        }
        fprintf(stderr, "OpenGL Debug Message:\n\tSeverity: %s\n\t"
                "Source: %s\n\tType: %s\n\tMessage: %s\n", severity_msg, src_msg,
                type_msg, message);
        assert(severity != GL_DEBUG_SEVERITY_HIGH && type != GL_DEBUG_TYPE_ERROR);
}


static void
mouse_button_cb(GLFWwindow *window, int button, int action, int mods)
{
        cur_quat   = mouse_quat(window);
        glfwGetCursorPos(window, &x_pos_prev, &y_pos_prev);
}


static void 
cursor_pos_cb(GLFWwindow *window, double x_pos, double y_pos)
{
        /* ImGui widget is hovered so we ingore all events */
        if (hover)
                return;

        /* left click does rotation of origin centered arcball */
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
                quat_t const tmp = mouse_quat(window);
                quat_t new_q;
                quat_mul(&new_q, &cur_quat, &tmp);
                cur_quat = tmp;

                mat4_t rotation_diff;
                quat_mat4(&rotation_diff, &new_q);
                mat4_mul(&rotation, &rotation_diff, &rotation);

                goto update;
        }

        /* right click does panning of the camera */
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
                int width, height;
                glfwGetWindowSize(window, &width, &height);

                /* less efficient but more clean */
                mat4_t panning_diff = MAT4;
                panning_diff.c[3][0] = 2.0f*(x_pos - x_pos_prev)/width;
                panning_diff.c[3][1] = -2.0f*(y_pos - y_pos_prev)/height;

                mat4_mul(&panning, &panning_diff, &panning);

                x_pos_prev = x_pos;
                y_pos_prev = y_pos;

                goto update;
        }

        return;
update:
        mat4_t tmp;
        mat4_mul(&tmp, &zoom, &rotation);
        mat4_mul(&ubo_mapped->view, &panning, &tmp);
}


static void 
scroll_cb(GLFWwindow *window, double xoffset, double yoffset)
{
        mat4_t zoom_diff = MAT4;
        if (yoffset > 0.0)
                mat4_scale(&zoom_diff, 1.1f, 1.1f, 1.1f);
        else
                mat4_scale(&zoom_diff, 0.9f, 0.9f, 0.9f);

        mat4_mul(&zoom, &zoom_diff, &zoom);

        mat4_t tmp;
        mat4_mul(&tmp, &zoom, &rotation);
        mat4_mul(&ubo_mapped->view, &panning, &tmp);
}

static void
key_cb(GLFWwindow *window, int key, int scancode, int action, int mods){
        if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE){
                glfwSetWindowShouldClose(window, 1);
                return;
        }

        ImGui_ImplGlfwGL3_KeyCallback(window, key, scancode, action, mods);
}

static void
realize()
{
        glDebugMessageCallback(debug, NULL);
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
        glEnable(GL_DEPTH_TEST);
        glPointSize(5.0f);

        int gl_ver_maj = 0, gl_ver_min = 0;
        glGetIntegerv(GL_MAJOR_VERSION, &gl_ver_maj);
        glGetIntegerv(GL_MINOR_VERSION, &gl_ver_min);
        printf("OpenGL Version %d.%d\n", gl_ver_maj, gl_ver_min);


        char *str = pk_io_read_file("data/vs.glsl");
        programs[0] = glCreateShaderProgramv(GL_VERTEX_SHADER, 1, &str);
        free(str);
        if (programs[0] == 0){
            char info_log[1024] = {0};
            GLsizei info_len = 0;
            glGetProgramInfoLog(programs[0], 1024, &info_len, info_log);
            printf("Program Failed to Compile or Link:\n%s\n--------", info_log);
        }

        str = pk_io_read_file("data/fs.glsl");
        programs[1] = glCreateShaderProgramv(GL_FRAGMENT_SHADER, 1, &str);
        free(str);
        if (programs[1] == 0){
            char info_log[1024] = {0};
            GLsizei info_len = 0;
            glGetProgramInfoLog(programs[0], 1024, &info_len, info_log);
            printf("Program Failed to Compile or Link:\n%s\n--------", info_log);
        }

        glCreateProgramPipelines(1, &pipeline);
        glBindProgramPipeline(pipeline);
        glUseProgramStages(pipeline, GL_VERTEX_SHADER_BIT, programs[0]);
        glUseProgramStages(pipeline, GL_FRAGMENT_SHADER_BIT, programs[1]);

        struct PdMeshSurface *mesh;
        if (mesh_filename){
                printf("Loading mesh %s\n", mesh_filename);
                char *str = pk_io_read_file(mesh_filename);
                assert(str);
                mesh = pd_mesh_surface_mk_from_json(str);
                free(str);
        } else
                mesh = pd_mesh_surface_mk_grid(resolution_x, resolution_y);

        assert(mesh);
        pd_mesh_print_info(mesh);

        triangles_count = mesh->n_indices;
        n_positions = mesh->n_positions;


        GLbitfield const flags = GL_MAP_COHERENT_BIT
                               | GL_MAP_PERSISTENT_BIT
                               | GL_MAP_WRITE_BIT;

        /* allocate buffers */
        glCreateBuffers(1, &vbo);
        glNamedBufferStorage(vbo, mesh->n_positions*3*sizeof *mesh->positions, mesh->positions, flags);

        /* map position buffer, topology of mesh does not change so we do not map indices */
        positions_mapped = (float *)glMapNamedBufferRange(vbo, 0, mesh->n_positions*3*sizeof *mesh->positions, flags);


        glCreateBuffers(N_EBOS, ebos);
        glNamedBufferStorage(ebos[EBO_TRIANGLES], mesh->n_indices*sizeof *mesh->indices, mesh->indices, 0);

        /* TODO: DOD on springs so we can pass pointer directly */
        lines_count             = 2*mesh->n_springs;
        size_t const lines_size = lines_count*sizeof *mesh->indices;
        uint32_t *lines_indices = (uint32_t *)malloc(lines_size);
        for (uint32_t i = 0; i < mesh->n_springs; ++i)
                memcpy(lines_indices + 2*i, mesh->springs[i].i, 2*sizeof *lines_indices);
        glNamedBufferStorage(ebos[EBO_LINES], lines_size, lines_indices, 0);
        free(lines_indices);


        glCreateBuffers(1, &ubo);
        glNamedBufferStorage(ubo, sizeof *ubo_mapped, NULL, flags | GL_MAP_READ_BIT);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo);

        ubo_mapped = (struct MatBlock *)glMapNamedBufferRange(ubo, 0, sizeof *ubo_mapped, flags | GL_MAP_READ_BIT);

        /* flip y and z axes */
        ubo_mapped->model         = MAT4;
        ubo_mapped->model.c[1][1] = 0.0f;
        ubo_mapped->model.c[1][2] = 1.0f;
        ubo_mapped->model.c[2][1] = 1.0f;
        ubo_mapped->model.c[2][2] = 0.0f;

        ubo_mapped->view       = MAT4;

        float const near = -5.0f;
        float const far  =  5.0f;
        ubo_mapped->projection = MAT4;
        ubo_mapped->projection.c[2][2] = -2.0f/(far - near);
        ubo_mapped->projection.c[3][2] = -(far + near)/(far - near);


        glCreateVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glVertexArrayAttribFormat(vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
        glVertexArrayAttribBinding(vao, 0, 0);
        glEnableVertexArrayAttrib(vao, 0);

        glBindVertexBuffer(0, vbo, 0, sizeof (GLfloat[3]));


        for (int i = 0; i < n_solvers; ++i) {
                mesh->attachments[0].position[0] = i;
                mesh->attachments[0].position[1] = 1.0f;
                mesh->attachments[0].position[2] = 0.0f;

                solvers[i] = pd_solver_alloc(mesh->positions, mesh->n_positions,
                                             mesh->attachments, mesh->n_attachments,
                                             mesh->springs, mesh->n_springs,
                                             timestep);
        }
        pd_mesh_surface_free(mesh);
}


static void
simulate()
{
        for (uint32_t i = 0; i < n_iterations; ++i)
                #pragma omp parallel for
                for (int j = 0; j < n_solvers; ++j)
                        pd_solver_advance(solvers[j]);

        /* reset the attachment constraint on shared border */
        float const *const pos[2] = {
                pd_solver_map_positions(solvers[0]),
                pd_solver_map_positions(solvers[1]),
        };

        struct PdConstraintAttachment *attachments[2] = {
                pd_solver_map_attachments(solvers[0]),
                pd_solver_map_attachments(solvers[1]),
        };

        for (uint32_t i = 1; i <= resolution_y; ++i)
                for (int j = 0; j < 3; ++j) {
                        attachments[0][i].position[j] = pos[1][3*(resolution_x*i - 2) + j];
                        attachments[1][i].position[j] = pos[0][3*(resolution_x*i - 2) + j];
                }
        
}


static void
render()
{
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        /* TODO: avoid memcpy by passing pointer to advance solver */
        for (int i = 0; i < n_solvers; ++i) {
                memcpy(positions_mapped, pd_solver_map_positions(solvers[i]), n_positions*3*sizeof *positions_mapped);

                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebos[EBO_TRIANGLES]);
                glDrawElements(GL_TRIANGLES, triangles_count, GL_UNSIGNED_INT, NULL);

                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebos[EBO_LINES]);
                glDrawElements(GL_LINES, lines_count, GL_UNSIGNED_INT, NULL);
                glDrawArrays(GL_POINTS, 0, n_positions);
                glFinish();
        }
}


static void
render_gui()
{
        ImGui_ImplGlfwGL3_NewFrame();

        ImGuiIO &io = ImGui::GetIO();
        ImGui::Text("Average %.3f ms/frame (%.1f FPS)", 1000.f / io.Framerate, io.Framerate);
        ImGui::Text("Solver: %s", pd_solver_name(solvers[0]));
        ImGui::Text("Local CMA/solve:  %.3fms", pd_solver_local_cma(solvers[0]));
        ImGui::Text("Global CMA/solve: %.3fms", pd_solver_global_cma(solvers[0]));

        if (ImGui::InputFloat3("gravity vector", gravity.v, -1, ImGuiInputTextFlags_CharsDecimal))
                for (int i = 0; i < n_solvers; ++i)
                        pd_solver_set_ext_force(solvers[i], gravity.v);

        pd_solver_draw_ui(solvers[0]);

        ImGui::Render();
}


static void
resize(GLFWwindow *window, int width, int height)
{
        glViewport(0, 0, width, height);

        float const aspect = (float)width/height;

        ubo_mapped->projection.c[0][0] = fmin(1.0f, 1.0f/aspect);
        ubo_mapped->projection.c[1][1] = fmin(1.0f, aspect);
}


static void
unrealize()
{
        for (int i = 0; i < n_solvers; ++i)
                pd_solver_free(solvers[i]);

        glDeleteProgram(programs[0]);
        glDeleteProgram(programs[1]);
        glDeleteProgramPipelines(1, &pipeline);
        glDeleteVertexArrays(1, &vao);
        glDeleteBuffers(2, ebos);
        glDeleteBuffers(1, &vbo);
        glDeleteBuffers(1, &ubo);
}

struct vec2 {
    int x, y;
};
std::istream& operator>>(std::istream &is, vec2 &v){
    is >> v.x >> v.y;
    return is;
}

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

int
main(int argc, char **argv)
{
        if (arg_flag(argv, argv + argc, "-h")){
            printf("Usage: ./pd_benchmark [options]\n"
                "\t--size <x> <y>       Cloth mesh size (default 10 10)\n"
                "\t--mesh <filename>    Tet mesh file to load\n"
                "\t-n <number>          Number of iterations of projective dynamics per timestep (default 10)\n");
            return 0;
        }
        if (arg_flag(argv, argv + argc, "--size")){
                std::array<int, 2> cloth_resolution = get_arg<int, 2>(argv, argv + argc, "--size");
                printf("Using cloth of size %dx%d\n", cloth_resolution[0], cloth_resolution[1]);
                resolution_x = cloth_resolution[0];
                resolution_y = cloth_resolution[1];
        } else if (arg_flag(argv, argv + argc, "--mesh")){
                mesh_filename = *(std::find(argv, argv + argc, std::string("--mesh")) + 1);
        }
        if (arg_flag(argv, argv + argc, "-n")){
                n_iterations = get_arg<uint32_t>(argv, argv + argc, "-n");
                if (n_iterations == 0){
                        printf("iteration count must be > 0! Forcing to 1\n");
                        n_iterations = 1;
                }
                timestep = 1.0f/(60.0f*n_iterations);
        }

        if (!glfwInit())
                exit(EXIT_FAILURE);

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

        GLFWwindow *window = glfwCreateWindow(800, 800, "Projective Dynamics",
                                              NULL, NULL);
        if (!window) {
                fputs("Failed to create window\n", stderr);
                glfwTerminate();
                exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);

        glfwSetCharCallback(window, ImGui_ImplGlfwGL3_CharCallback);
        glfwSetCursorPosCallback(window, cursor_pos_cb);
        glfwSetKeyCallback(window, key_cb);
        glfwSetMouseButtonCallback(window, mouse_button_cb);
        glfwSetScrollCallback(window, scroll_cb);
        glfwSetWindowSizeCallback(window, resize);

        ImGui_ImplGlfwGL3_Init(window, false);

        realize();

        while (!glfwWindowShouldClose(window)) {
                simulate();
                render();
                render_gui();

                hover = ImGui::IsMouseHoveringAnyWindow();

                glfwSwapBuffers(window);
                glfwPollEvents();
        }

        unrealize();

        glfwDestroyWindow(window);
        glfwTerminate();
}
