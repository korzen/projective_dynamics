#ifndef PD_MESH_H
#define PD_MESH_H


#include <assert.h>
#include <tgmath.h>

#include <jsmn/jsmn.h>

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
        struct PdConstraintSpring s = { i0, i1, rest_length, };
        return s;
}


/* grid generated in NDC coordinates on xz-plane with corner at (1, 0, 1) */
struct PdMeshSurface *
pd_mesh_surface_mk_grid(uint32_t const n_x, uint32_t const n_y)
{
        /* must generate at least one quad (n_x = 2, n_y = 2) */
        assert(n_x > 1 && n_y > 1);

        struct PdMeshSurface *m = (struct PdMeshSurface *)malloc(sizeof *m);

        m->n_positions = n_x*n_y;
        m->positions   = (float *)malloc(3*m->n_positions*sizeof *m->positions);
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
        m->indices   = (uint32_t *)malloc(m->n_indices*sizeof *m->indices);
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
        m->springs   = (struct PdConstraintSpring *)malloc(n_springs*sizeof *m->springs);
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

        /* create net */
        m->n_attachments  = 2;
        m->attachments    = (struct PdConstraintAttachment *)malloc(m->n_attachments*sizeof *m->attachments);

        /*
        m->attachments[0] = (struct PdConstraintAttachment){ 0, { m->positions[0], m->positions[1], m->positions[2], }, };
        m->attachments[1] = (struct PdConstraintAttachment){ (n_x - 1), { m->positions[3*(n_x - 1)], m->positions[3*(n_x - 1) + 1], m->positions[3*(n_x - 1) + 2], }, };
        */
        m->attachments[0] = (struct PdConstraintAttachment){ 0, { -1.0f, -1.0f, 1.0f, }, };
        m->attachments[1] = (struct PdConstraintAttachment){ (n_x - 5), { 1.0f, -1.0f, 1.0f, }, };

        //m->attachments[2] = (struct PdConstraintAttachment){ (n_y - 1)*n_x, { -1.0f, 1.0f, 1.0f, }, };
        //m->attachments[3] = (struct PdConstraintAttachment){ (n_x - 1) + (n_y - 1)*n_x, { 1.0f, 1.0f, 1.0f, }, };


        return m;
}


static int
json_eq(char const *json, jsmntok_t const *tok, char const *str)
{
        return tok->type == JSMN_STRING && (size_t)tok->end - (size_t)tok->start == strlen(str) && !strncmp(json + tok->start, str, tok->end - tok->start);
}


static int
parse_attachment(struct PdConstraintAttachment *a, char const *json, jsmntok_t *t, int i, int const r)
{
        for (int done = 0; i < r && done < 2; ++i) {
                if (json_eq(json, t + i, "i")) {
                        a->i = strtoul(json + t[++i].start, NULL, 10);
                        ++done;
                } else if (json_eq(json, t + i, "position")) {
                        ++i;
                        for (int j = 0; j < 3; ++j)
                                a->position[j] = strtof(json + t[++i].start, NULL);
                        ++done;
                }
        }
        return i;
}


static int
parse_spring(struct PdConstraintSpring *a, char const *json, jsmntok_t *t, int i, int const r)
{
        for (int done = 0; i < r && done < 2; ++i) {
                if (json_eq(json, t + i, "i")) {
                        ++i;
                        for (int j = 0; j < 2; ++j)
                                a->i[j] = strtoul(json + t[++i].start, NULL, 10);
                        ++done;
                } else if (json_eq(json, t + i, "rest_length")) {
                        a->rest_length = strtof(json + t[++i].start, NULL);
                        ++done;
                }
        }
        return i;
}


struct PdMeshSurface *
pd_mesh_surface_mk_from_json(char const *json)
{
        jsmn_parser parser;
        jsmn_init(&parser);

        size_t const n_tokens = 1 << 24;
        jsmntok_t *t = (jsmntok_t *)malloc(n_tokens*sizeof *t);
        int const r = jsmn_parse(&parser, json, strlen(json), t, n_tokens);
        if (r < 0) {
                puts("problem parsing the data");
                puts("not enough tokens :/");
                return NULL;
        }

        struct PdMeshSurface *m = (struct PdMeshSurface *)malloc(sizeof *m);
        m->n_positions   = 0;
        m->n_indices     = 0;
        m->n_attachments = 0;
        m->n_springs     = 0;


        for (int i = 0; i < r;) {
                if (json_eq(json, t + i, "attachments")) {
                        m->n_attachments = t[++i].size;
                        m->attachments   = (struct PdConstraintAttachment *)malloc(m->n_attachments*sizeof *m->attachments);

                        for (uint32_t o = 0; o < m->n_attachments; ++o)
                                i = parse_attachment(m->attachments + o, json, t, i, r);
                } else if (json_eq(json, t + i, "indices")) {
                        m->n_indices = 3*t[++i].size;
                        m->indices   = (uint32_t *)malloc(m->n_indices*sizeof *m->indices);

                        uint32_t o = 0;
                        for (++i; i < r && t[i].type == JSMN_ARRAY; ++i) {
                                for (int j = 0; j < 3; ++j, ++o)
                                        m->indices[o] = strtoul(json + t[++i].start, NULL, 10);
                        }
                } else if (json_eq(json, t + i, "springs")) {
                        m->n_springs = t[++i].size;
                        m->springs   = (struct PdConstraintSpring *)malloc(m->n_springs*sizeof *m->springs);

                        for (uint32_t o = 0; o < m->n_springs; ++o)
                                i = parse_spring(m->springs + o, json, t, i, r);
                } else if (json_eq(json, t + i, "vertices")) {
                        m->n_positions = t[++i].size;
                        m->positions   = (float *)malloc(3*m->n_positions*sizeof *m->positions);

                        uint32_t o = 0;
                        for (++i; i < r && t[i].type == JSMN_ARRAY; ++i) {
                                for (int j = 0; j < 3; ++j, ++o)
                                        m->positions[o] = strtof(json + t[++i].start, NULL);
                        }
                } else
                        ++i;
        }
        free(t);

        return m;
}

/* Binary mesh file format:
 * uint32_t: # of attachments
 * uint32_t: # of indices
 * uint32_t: # of springs
 * uint32_t: # of vertices
 * [attachments...]: each attachment should be:
 *      int32_t: index of the attachment
 *      [f32, f32, f32]: position that is attached too
 * [indices...]: indices for each triangle, each should be:
 *      [int32_t, int32_t, int32_t]: first, second, third vert for each tri
 * [springs...]: spring attachments, each should be:
 *      [int32_t, int32_t]: indices of the attached vertices
 *      f32: rest length of the spring
 * [vertices...]: vertex positions, each should be:
 *      [f32, f32, f32]: x, y, z coordinates of the vertex
 */
struct PdMeshSurface *
pd_mesh_surface_mk_from_binary(char const *fname){
        FILE *fp = fopen(fname, "r");
        if (!fp){
                printf("Failed to open file %s\n", fname);
                return NULL;
        }
        struct PdMeshSurface *m = (struct PdMeshSurface *)malloc(sizeof *m);
        m->n_attachments = 0;
        m->n_indices     = 0;
        m->n_springs     = 0;
        m->n_positions   = 0;
        // Read counts of the various mesh components
        if (fread(&m->n_attachments, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to read # of attachments\n");
                free(m);
                fclose(fp);
                return NULL;
        }
        printf("read n_attachments = %u\n", m->n_attachments);
        if (fread(&m->n_indices, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to read # of indices\n");
                free(m);
                fclose(fp);
                return NULL;
        }
        if (fread(&m->n_springs, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to read # of springs\n");
                free(m);
                fclose(fp);
                return NULL;
        }
        if (fread(&m->n_positions, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to read # of vertices\n");
                free(m);
                fclose(fp);
                return NULL;
        }
        // Read the attachments
        m->attachments = (struct PdConstraintAttachment *)malloc(m->n_attachments*sizeof *m->attachments);
        if (fread(m->attachments, sizeof *m->attachments, m->n_attachments, fp) != m->n_attachments){
                printf("Failed to read attachments\n");
                free(m->attachments);
                free(m);
                fclose(fp);
                return NULL;
        }
        // Read the indices
        m->indices = (uint32_t *)malloc(m->n_indices*sizeof *m->indices);
        if (fread(m->indices, sizeof *m->indices, m->n_indices, fp) != m->n_indices){
                printf("Failed to read indices\n");
                free(m->indices);
                free(m->attachments);
                free(m);
                fclose(fp);
                return NULL;
        }
        // Read the springs
        m->springs = (struct PdConstraintSpring *)malloc(m->n_springs*sizeof *m->springs);
        if (fread(m->springs, sizeof *m->springs, m->n_springs, fp) != m->n_springs){
                printf("Failed to read springs\n");
                free(m->springs);
                free(m->indices);
                free(m->attachments);
                free(m);
                fclose(fp);
                return NULL;
        }
        // Read the vertices
        m->positions = (float *)malloc(3*m->n_positions*sizeof *m->positions);
        if (fread(m->positions, sizeof *m->positions, 3 * m->n_positions, fp) != 3 * m->n_positions){
                printf("Failed to read positions\n");
                free(m->positions);
                free(m->springs);
                free(m->indices);
                free(m->attachments);
                free(m);
                fclose(fp);
                return NULL;
        }
        fclose(fp);
        return m;
}

// Write the mesh to the file in our binary format, returns 1 if success and 0 if failure
int
pd_mesh_surface_write_binary(struct PdMeshSurface const *m, char const *fname){
        FILE *fp = fopen(fname, "w");
        if (!fp){
                printf("Failed to open file %s\n", fname);
                return 0;
        }
        // Write counts of the various mesh components
        if (fwrite(&m->n_attachments, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to write # of attachments\n");
                fclose(fp);
                return 0;
        }
        if (fwrite(&m->n_indices, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to write # of indices\n");
                fclose(fp);
                return 0;
        }
        if (fwrite(&m->n_springs, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to write # of springs\n");
                fclose(fp);
                return 0;
        }
        if (fwrite(&m->n_positions, sizeof(uint32_t), 1, fp) != 1){
                printf("Failed to write # of vertices\n");
                fclose(fp);
                return 0;
        }
        // write the attachments
        if (fwrite(m->attachments, sizeof *m->attachments, m->n_attachments, fp) != m->n_attachments){
                printf("Failed to write attachments\n");
                fclose(fp);
                return 0;
        }
        // write the indices, note: we have 3 ints per index
        if (fwrite(m->indices, sizeof *m->indices, m->n_indices, fp) != m->n_indices){
                printf("Failed to write indices\n");
                fclose(fp);
                return 0;
        }
        // write the springs
        if (fwrite(m->springs, sizeof *m->springs, m->n_springs, fp) != m->n_springs){
                printf("Failed to write springs\n");
                fclose(fp);
                return 0;
        }
        // write the vertices
        if (fwrite(m->positions, sizeof *m->positions, 3 * m->n_positions, fp) != 3 * m->n_positions){
                printf("Failed to write positions\n");
                fclose(fp);
                return 0;
        }
        fclose(fp);
        return 1;
}

void
pd_mesh_print_info(struct PdMeshSurface const *m)
{
        printf("n_attachments %u\n", m->n_attachments);
        printf("n_indices     %u\n", m->n_indices);
        printf("n_positions   %u\n", m->n_positions);
        printf("n_springs     %u\n", m->n_springs);
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
