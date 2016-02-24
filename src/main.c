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

#include "pd_mesh.h"
#include "pd_solver.h"


static GLsizei count;
static GLuint  pipeline;
static GLuint  programs[2];
static GLuint  vao;
static GLuint  vbos[2];

static uint32_t         n_positions;
static float           *positions_mapped;
static struct PdSolver *solver;


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
        gdk_gl_context_set_required_version(context, 3, 3);

        return context;
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


        struct PdMeshSurface *mesh = pd_mesh_surface_mk_grid(16, 8);
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


        /* set up GL window */
        GtkWidget *gl_area = gtk_gl_area_new();
        gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(gl_area), TRUE);

        g_signal_connect(gl_area, "create-context",
                         G_CALLBACK(create_context), NULL);
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
