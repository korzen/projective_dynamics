#ifndef PD_MESH_H
#define PD_MESH_H


#include <assert.h>
#include <tgmath.h>

#include "pd_constraint.h"


struct PdMeshSurface {
        float *positions;
        uint32_t n_positions;

        uint32_t *indices; /* we want to go above 65k */
        uint32_t n_indices;

        struct PdConstraintAttachment *attachments;
        uint32_t n_attachments;

        struct PdConstraintSpring *springs;
        uint32_t n_springs;

};


static struct PdConstraintSpring
mk_spring(float const *positions, uint32_t const i0, uint32_t const i1)
{
        float const *p0 = positions + 3*i0;
        float const *p1 = positions + 3*i1;
        float const rest_length = sqrt(pow(p1[0] - p0[0], 2) + pow(p1[1] - p0[1], 2) + pow(p1[2] - p0[2], 2));
        struct PdConstraintSpring s = {
                .i[0] = i0,
                .i[1] = i1,
                .rest_length = rest_length,
        };
        return s;
}


/* grid generated in NDC coordinates on xz-plane with corner at (1, 0, 1) */
struct PdMeshSurface *
pd_mesh_surface_mk_grid(uint32_t const n_x, uint32_t const n_y)
{
        /* must generate at least one quad (n_x = 2, n_y = 2) */
        assert(n_x > 1 && n_y > 1);

        struct PdMeshSurface *m = malloc(sizeof *m);

        m->n_positions = n_x*n_y;
        m->positions   = malloc(3*m->n_positions*sizeof *m->positions);
        for (uint32_t j = 0; j < n_y; ++j) {
                for (uint32_t i = 0; i < n_x; ++i) {
                        /* convert to NDC */
                        float *p = m->positions + 3*(i + j*n_x);
                        p[0] = 2.0f*i/(n_x - 1.0f) - 1.0f;
                        p[1] = 0.0f;
                        p[2] = 1.0f - 2.0f*j/(n_y - 1.0f);
                }
        }

        /* diagonal goes from top-left corner */
        /* 3 indices per triangles; 2 triangles per quad */
        m->n_indices = 3*2*(n_x - 1)*(n_y - 1);
        m->indices   = malloc(m->n_indices*sizeof *m->indices);
        for (uint32_t j = 0; j + 1 < n_y; ++j) {
                for (uint32_t i = 0; i + 1 < n_x; ++i) {
                        /* generate two triangles for the quad */
                        uint32_t const quad[4] = {
                                i + j*n_x,
                                i + (j + 1)*n_x,
                                (i + 1) + (j + 1)*n_x,
                                (i + 1) + j*n_x,
                        };

                        /* TODO: use loop and tessellate in it */
                        uint32_t *t = m->indices + 3*2*(i + j*(n_x - 1));
                        t[0] = quad[0];
                        t[1] = quad[1];
                        t[2] = quad[2];

                        t += 3;
                        t[0] = quad[0];
                        t[1] = quad[2];
                        t[2] = quad[3];
                }
        }

        /* spring constraints of quad mesh */
        uint32_t const n_springs = 2*(n_x - 1)*(n_y - 1) + (n_y - 1) + (n_x - 1);
        m->n_springs = 0;
        m->springs   = malloc(n_springs*sizeof *m->springs);
        /* generate upper corner constraints */
        for (uint32_t j = 0; j + 1 < n_y; ++j) {
                for (uint32_t i = 0; i + 1 < n_x; ++i) {
                        m->springs[m->n_springs++] = mk_spring(m->positions, i + j*n_x, (i + 1) + j*n_x);
                        m->springs[m->n_springs++] = mk_spring(m->positions, i + j*n_x, i + (j + 1)*n_x);
                }
        }
        /* generate right edge constraints */
        for (uint32_t j = 0; j + 1 < n_y; ++j)
                m->springs[m->n_springs++] = mk_spring(m->positions, j*n_x + (n_x - 1), (j + 1)*n_x + (n_x - 1));

        /* generate bottom edge constraints */
        for (uint32_t i = 0; i + 1 < n_x; ++i)
                m->springs[m->n_springs++] = mk_spring(m->positions, (n_y - 1)*n_x + i, (n_y - 1)*n_x + (i + 1));

        m->attachments   = NULL;
        m->n_attachments = 0;

        return m;
}


void
pd_mesh_surface_free(struct PdMeshSurface *m)
{
        free(m->positions);
        free(m->indices);
        free(m->attachments);
        free(m->springs);
        free(m);
}


#endif
