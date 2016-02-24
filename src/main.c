#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>

#include <epoxy/gl.h>
#include <gtk/gtk.h>
#include <pk/pk_io.h>
#include <pk/pk_linalg.h>

#include "pd_mesh.h"
#include "pd_solver.h"


static GLsizei count;
static GLuint  pipeline;
static GLuint  programs[2];
static GLuint  ubo;
static GLuint  vao;
static GLuint  vbos[2];

static quat_t cur_quat;
static uint32_t n_positions;
static float *positions_mapped;
static struct {
        mat4_t model;
        mat4_t view;
        mat4_t projection;
} *ubo_mapped;
static struct PdSolver *solver;


static quat_t
mouse_quat(GtkWidget *widget, GdkEvent *event)
{
        gdouble x_win, y_win;
        gdk_event_get_coords(event, &x_win, &y_win);

        int const width  = gtk_widget_get_allocated_width(widget);
        int const height = gtk_widget_get_allocated_height(widget);

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
        fprintf(stderr, "%s\n", message);
}


static GdkGLContext *
create_context(GtkGLArea *area)
{
        GdkWindow *window = gtk_widget_get_window(GTK_WIDGET(area));

        GError *error = NULL;
        GdkGLContext *context = gdk_window_create_gl_context(window, &error);
        if (error) {
                fprintf(stderr, "%s\n", error->message);
                g_error_free(error);

                return NULL;
        }

        gdk_gl_context_set_debug_enabled(context, TRUE);
        gdk_gl_context_set_forward_compatible(context, TRUE);
        gdk_gl_context_set_required_version(context, 4, 5);

        return context;
}


static gboolean
esc_key_press_event(GtkWidget *widget, GdkEvent *event, gpointer user_data)
{
        guint keyval;
        gdk_event_get_keyval(event, &keyval);

        if (keyval == GDK_KEY_Escape)
                gtk_main_quit();

        return TRUE;
}


static gboolean
key_press_event(GtkWidget *widget, GdkEvent *event, gpointer user_data)
{
        guint keyval;
        gdk_event_get_keyval(event, &keyval);

        switch (keyval) {
        case GDK_KEY_Left:
                positions_mapped[0] -= 0.01f;
                goto redraw;
        case GDK_KEY_Right:
                positions_mapped[0] += 0.01;
                goto redraw;
        default:
                return FALSE;
        }

redraw:
        gtk_gl_area_queue_render(user_data);
        return TRUE;
}


static gboolean
button_press_event(GtkWidget *widget, GdkEvent *event, gpointer user_data)
{
        cur_quat = mouse_quat(widget, event);

        return TRUE;
}


static gboolean
motion_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data)
{
        GdkModifierType state;
        gdk_event_get_state(event, &state);

        if (!(state & GDK_BUTTON1_MASK))
                return FALSE;

        quat_t const tmp = mouse_quat(widget, event);
        quat_t new;
        quat_mul(&new, &cur_quat, &tmp);
        cur_quat = tmp;

        mat4_t rotation;
        quat_mat4(&rotation, &new);
        mat4_mul(&ubo_mapped->view, &rotation, &ubo_mapped->view);

        gtk_gl_area_queue_render(GTK_GL_AREA(widget));

        return TRUE;
}


static void
realize(GtkWidget *widget, gpointer user_data)
{
        gtk_gl_area_make_current(GTK_GL_AREA(widget));

        glDebugMessageCallback(debug, NULL);

        char *str = pk_io_read_file("data/vs.glsl");
        programs[0] = glCreateShaderProgramv(GL_VERTEX_SHADER, 1, &str);
        free(str);

        str = pk_io_read_file("data/fs.glsl");
        programs[1] = glCreateShaderProgramv(GL_FRAGMENT_SHADER, 1, &str);
        free(str);

        glCreateProgramPipelines(1, &pipeline);
        glUseProgramStages(pipeline, GL_VERTEX_SHADER_BIT, programs[0]);
        glUseProgramStages(pipeline, GL_FRAGMENT_SHADER_BIT, programs[1]);
        glBindProgramPipeline(pipeline);


        struct PdMeshSurface *mesh = pd_mesh_surface_mk_grid(16, 16);
        count = mesh->n_indices;
        n_positions = mesh->n_positions;


        GLbitfield const flags = GL_MAP_COHERENT_BIT
                               | GL_MAP_PERSISTENT_BIT
                               | GL_MAP_WRITE_BIT;

        /* allocate buffers */
        glCreateBuffers(2, vbos);
        glNamedBufferStorage(vbos[0], mesh->n_positions*3*sizeof *mesh->positions, mesh->positions, flags);
        glNamedBufferStorage(vbos[1], mesh->n_indices*sizeof *mesh->indices, mesh->indices, 0);

        /* map position buffer, topology of mesh does not change */
        positions_mapped = glMapNamedBufferRange(vbos[0], 0, mesh->n_positions*sizeof *mesh->positions, flags);


        glCreateBuffers(1, &ubo);
        glNamedBufferStorage(ubo, sizeof *ubo_mapped, NULL, flags | GL_MAP_READ_BIT);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo);

        ubo_mapped = glMapNamedBufferRange(ubo, 0, sizeof *ubo_mapped, flags | GL_MAP_READ_BIT);

        /* flip y and z axes */
        ubo_mapped->model         = MAT4;
        ubo_mapped->model.c[1][1] = 0.0f;
        ubo_mapped->model.c[1][2] = 1.0f;
        ubo_mapped->model.c[2][1] = 1.0f;
        ubo_mapped->model.c[2][2] = 0.0f;

        ubo_mapped->view       = MAT4;
        ubo_mapped->projection = MAT4;
        ubo_mapped->projection.c[2][2] = -1.0f;


        glCreateVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glVertexArrayAttribFormat(vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
        glVertexArrayAttribBinding(vao, 0, 0);
        glEnableVertexArrayAttrib(vao, 0);

        glBindVertexBuffer(0, vbos[0], 0, sizeof (GLfloat[3]));
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos[1]);


        solver = pd_solver_alloc(mesh->positions, mesh->n_positions,
                                 mesh->attachments, mesh->n_attachments,
                                 mesh->springs, mesh->n_springs);
        pd_mesh_surface_free(mesh);
}


static gboolean
render(GtkGLArea *area, GdkGLContext *context, gpointer user_data)
{
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        /* TODO: avoid memcpy by passing pointer to advance solver */
        float const timestep = 1.0f/60.0f;
        pd_solver_advance(solver, timestep);

        memcpy(positions_mapped, pd_solver_map_positions(solver), n_positions*3*sizeof *positions_mapped);

        /*
        glDrawElements(GL_TRIANGLES, (sizeof indices)/(sizeof indices[0]),
                       GL_UNSIGNED_BYTE, NULL);
        */
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glPointSize(5.0f);
        glDrawArrays(GL_POINTS, 0, n_positions);
        glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT, NULL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        return TRUE;
}


static void
unrealize(GtkWidget *widget, gpointer user_data)
{
        printf("destructing\n");
        gtk_gl_area_make_current(GTK_GL_AREA(widget));

        glDeleteProgram(programs[0]);
        glDeleteProgram(programs[1]);
        glDeleteProgramPipelines(1, &pipeline);
        glDeleteVertexArrays(1, &vao);
}


gboolean
animate(GtkWidget *widget, GdkFrameClock *frame_clock, gpointer user_data)
{
        gtk_gl_area_queue_render(GTK_GL_AREA(widget));
        return G_SOURCE_CONTINUE;
}


int
main(int argc, char **argv)
{
        gtk_init(&argc, &argv);

        GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
        gtk_widget_add_events(window, GDK_KEY_PRESS_MASK);
        g_signal_connect(window, "key-press-event", G_CALLBACK(esc_key_press_event), NULL);

        /* set up GL window */
        GtkWidget *gl_area = gtk_gl_area_new();
        gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(gl_area), TRUE);
        gtk_widget_add_events(gl_area,
                              GDK_BUTTON_PRESS_MASK |
                              GDK_BUTTON_RELEASE_MASK |
                              GDK_POINTER_MOTION_MASK);

        g_signal_connect(gl_area, "button-press-event", G_CALLBACK(button_press_event), NULL);
        g_signal_connect(gl_area, "create-context", G_CALLBACK(create_context), NULL);
        g_signal_connect(gl_area, "motion-notify-event", G_CALLBACK(motion_notify_event), NULL);
        g_signal_connect(gl_area, "realize", G_CALLBACK(realize), NULL);
        g_signal_connect(gl_area, "render", G_CALLBACK(render), NULL);
        g_signal_connect(gl_area, "unrealize", G_CALLBACK(unrealize), NULL);

        gtk_container_add(GTK_CONTAINER(window), gl_area);
        gtk_widget_add_tick_callback(gl_area, animate, NULL, NULL);

        g_signal_connect(window, "key-press-event",
                         G_CALLBACK(key_press_event), gl_area);

        gtk_widget_show_all(window);

        gtk_main();
}
