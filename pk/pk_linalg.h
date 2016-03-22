#ifndef PK_LINALG_H
#define PK_LINALG_H

typedef union {
        struct { float x, y; };
        float v[2];
} vec2_t;

typedef union {
        struct { float x, y, z; };
        struct { vec2_t xy; float z_; };
        float v[3];
} vec3_t;

typedef union {
        struct { float x, y, z, w; };
        struct { vec3_t xyz; float w_; };
        float v[4];
} vec4_t;

/* column-major storage layout (glsl) */
typedef union {
        struct { vec4_t v0, v1, v2, v3; };
        float v[16];
        float c[4][4];
} mat4_t;

typedef union {
        struct { float w, x, y, z; };
        struct { float w_; vec3_t xyz; };
        float q[4];
} quat_t;


#define VEC2(x, y) (vec2_t){{x, y}}

#define VEC4(x, y, z, w) (vec4_t){{x, y, z, w}}

#define MAT4 (mat4_t){{ \
        {{ 1.0f, 0.0f, 0.0f, 0.0f, }}, \
        {{ 0.0f, 1.0f, 0.0f, 0.0f, }}, \
        {{ 0.0f, 0.0f, 1.0f, 0.0f, }}, \
        {{ 0.0f, 0.0f, 0.0f, 1.0f, }}, \
}}

#define QUAT(w, x, y, z) (quat_t){{w, x, y, z}}


vec2_t
vec2_add(vec2_t const u, vec2_t const v);

vec2_t
vec2_sub(vec2_t const u, vec2_t const v);

float
vec2_dot(vec2_t const u, vec2_t const v);

vec2_t
vec2_mul(float const s, vec2_t const v);

void
mat4_inv(mat4_t *restrict res, mat4_t const *restrict m);

void
mat4_mul(mat4_t *res, mat4_t const *m0, mat4_t const *m1);

void
mat4_scale(mat4_t *m, float const x, float const y, float const z);

void
quat_mul(quat_t *res, quat_t const *restrict q0, quat_t const *restrict q1);

void
quat_mat4(mat4_t *restrict m, quat_t const *restrict q);


#ifdef PK_LINALG_IMPLEMENTATION

#include <assert.h>
#include <tgmath.h>


vec2_t
vec2_add(vec2_t const u, vec2_t const v)
{
        return VEC2(u.x + v.x, u.y + v.y);
}


float
vec2_dot(vec2_t const u, vec2_t const v)
{
        return u.x*v.x + u.y*v.y;
}


float
vec2_length(vec2_t const v)
{
        return sqrtf(v.x*v.x + v.y*v.y);
}


vec2_t
vec2_mul(float const s, vec2_t const v)
{
        return VEC2(s*v.x, s*v.y);
}


vec2_t
vec2_normalize(vec2_t const v)
{
        return vec2_mul(1.0f/vec2_length(v), v);
}


vec2_t
vec2_sub(vec2_t const u, vec2_t const v)
{
        return VEC2(u.x - v.x, u.y - v.y);
}


void
vec4_add(vec4_t *restrict res, vec4_t const *restrict u, vec4_t const *restrict v)
{
        res->x = u->x + v->x;
        res->y = u->y + v->y;
        res->z = u->z + v->z;
        res->w = u->w + v->w;
}


void
vec4_sub(vec4_t *restrict res, vec4_t const *restrict u, vec4_t const *restrict v)
{
        res->x = u->x - v->x;
        res->y = u->y - v->y;
        res->z = u->z - v->z;
        res->w = u->w - v->w;
}


void
vec4_cross(vec4_t *restrict res, vec4_t const *restrict u, vec4_t const *restrict v)
{
        res->x = u->y*v->z - u->z*v->y;
        res->y = u->z*v->x - u->x*v->z;
        res->z = u->x*v->y - u->y*v->x;
        res->w = 0.0f;
}


float
vec4_dot(vec4_t const *restrict u, vec4_t const *restrict v)
{
        return u->x*v->x + u->y*v->y + u->z*v->z + u->w*v->w;
}


float
vec4_length(vec4_t const *v)
{
        return sqrtf(vec4_dot(v, v));
}


void
vec4_mul(vec4_t *restrict res, float const s, vec4_t const *restrict v)
{
        res->x = s*v->x;
        res->y = s*v->y;
        res->z = s*v->z;
        res->w = s*v->w;
}


void
vec4_normalize(vec4_t *v)
{
        vec4_mul(v, 1.0f/vec4_length(v), v);
}


/* lifted from MESA */
void
mat4_inv(mat4_t *restrict res, mat4_t const *restrict m)
{
    res->v[0] = m->v[5]  * m->v[10] * m->v[15] -
             m->v[5]  * m->v[11] * m->v[14] -
             m->v[9]  * m->v[6]  * m->v[15] +
             m->v[9]  * m->v[7]  * m->v[14] +
             m->v[13] * m->v[6]  * m->v[11] -
             m->v[13] * m->v[7]  * m->v[10];

    res->v[4] = -m->v[4]  * m->v[10] * m->v[15] +
              m->v[4]  * m->v[11] * m->v[14] +
              m->v[8]  * m->v[6]  * m->v[15] -
              m->v[8]  * m->v[7]  * m->v[14] -
              m->v[12] * m->v[6]  * m->v[11] +
              m->v[12] * m->v[7]  * m->v[10];

    res->v[8] = m->v[4]  * m->v[9] * m->v[15] -
             m->v[4]  * m->v[11] * m->v[13] -
             m->v[8]  * m->v[5] * m->v[15] +
             m->v[8]  * m->v[7] * m->v[13] +
             m->v[12] * m->v[5] * m->v[11] -
             m->v[12] * m->v[7] * m->v[9];

    res->v[12] = -m->v[4]  * m->v[9] * m->v[14] +
               m->v[4]  * m->v[10] * m->v[13] +
               m->v[8]  * m->v[5] * m->v[14] -
               m->v[8]  * m->v[6] * m->v[13] -
               m->v[12] * m->v[5] * m->v[10] +
               m->v[12] * m->v[6] * m->v[9];

    res->v[1] = -m->v[1]  * m->v[10] * m->v[15] +
              m->v[1]  * m->v[11] * m->v[14] +
              m->v[9]  * m->v[2] * m->v[15] -
              m->v[9]  * m->v[3] * m->v[14] -
              m->v[13] * m->v[2] * m->v[11] +
              m->v[13] * m->v[3] * m->v[10];

    res->v[5] = m->v[0]  * m->v[10] * m->v[15] -
             m->v[0]  * m->v[11] * m->v[14] -
             m->v[8]  * m->v[2] * m->v[15] +
             m->v[8]  * m->v[3] * m->v[14] +
             m->v[12] * m->v[2] * m->v[11] -
             m->v[12] * m->v[3] * m->v[10];

    res->v[9] = -m->v[0]  * m->v[9] * m->v[15] +
              m->v[0]  * m->v[11] * m->v[13] +
              m->v[8]  * m->v[1] * m->v[15] -
              m->v[8]  * m->v[3] * m->v[13] -
              m->v[12] * m->v[1] * m->v[11] +
              m->v[12] * m->v[3] * m->v[9];

    res->v[13] = m->v[0]  * m->v[9] * m->v[14] -
              m->v[0]  * m->v[10] * m->v[13] -
              m->v[8]  * m->v[1] * m->v[14] +
              m->v[8]  * m->v[2] * m->v[13] +
              m->v[12] * m->v[1] * m->v[10] -
              m->v[12] * m->v[2] * m->v[9];

    res->v[2] = m->v[1]  * m->v[6] * m->v[15] -
             m->v[1]  * m->v[7] * m->v[14] -
             m->v[5]  * m->v[2] * m->v[15] +
             m->v[5]  * m->v[3] * m->v[14] +
             m->v[13] * m->v[2] * m->v[7] -
             m->v[13] * m->v[3] * m->v[6];

    res->v[6] = -m->v[0]  * m->v[6] * m->v[15] +
              m->v[0]  * m->v[7] * m->v[14] +
              m->v[4]  * m->v[2] * m->v[15] -
              m->v[4]  * m->v[3] * m->v[14] -
              m->v[12] * m->v[2] * m->v[7] +
              m->v[12] * m->v[3] * m->v[6];

    res->v[10] = m->v[0]  * m->v[5] * m->v[15] -
              m->v[0]  * m->v[7] * m->v[13] -
              m->v[4]  * m->v[1] * m->v[15] +
              m->v[4]  * m->v[3] * m->v[13] +
              m->v[12] * m->v[1] * m->v[7] -
              m->v[12] * m->v[3] * m->v[5];

    res->v[14] = -m->v[0]  * m->v[5] * m->v[14] +
               m->v[0]  * m->v[6] * m->v[13] +
               m->v[4]  * m->v[1] * m->v[14] -
               m->v[4]  * m->v[2] * m->v[13] -
               m->v[12] * m->v[1] * m->v[6] +
               m->v[12] * m->v[2] * m->v[5];

    res->v[3] = -m->v[1] * m->v[6] * m->v[11] +
              m->v[1] * m->v[7] * m->v[10] +
              m->v[5] * m->v[2] * m->v[11] -
              m->v[5] * m->v[3] * m->v[10] -
              m->v[9] * m->v[2] * m->v[7] +
              m->v[9] * m->v[3] * m->v[6];

    res->v[7] = m->v[0] * m->v[6] * m->v[11] -
             m->v[0] * m->v[7] * m->v[10] -
             m->v[4] * m->v[2] * m->v[11] +
             m->v[4] * m->v[3] * m->v[10] +
             m->v[8] * m->v[2] * m->v[7] -
             m->v[8] * m->v[3] * m->v[6];

    res->v[11] = -m->v[0] * m->v[5] * m->v[11] +
               m->v[0] * m->v[7] * m->v[9] +
               m->v[4] * m->v[1] * m->v[11] -
               m->v[4] * m->v[3] * m->v[9] -
               m->v[8] * m->v[1] * m->v[7] +
               m->v[8] * m->v[3] * m->v[5];

    res->v[15] = m->v[0] * m->v[5] * m->v[10] -
              m->v[0] * m->v[6] * m->v[9] -
              m->v[4] * m->v[1] * m->v[10] +
              m->v[4] * m->v[2] * m->v[9] +
              m->v[8] * m->v[1] * m->v[6] -
              m->v[8] * m->v[2] * m->v[5];

    float const det = m->v[0]*res->v[0] + m->v[1]*res->v[4] + m->v[2]*res->v[8] + m->v[3]*res->v[12];

    assert(det != 0.0f);
    float const det_inv = 1.0f/det;

    for (int i = 0; i < 16; i++)
        res->v[i] = res->v[i]*det_inv;
}


void
mat4_scale(mat4_t *m, float const x, float const y, float const z)
{
        m->v[0] = x;
        m->v[5] = y;
        m->v[10] = z;
}


void
mat4_translate(mat4_t *m, float const x, float const y, float const z)
{
        m->v[12] = x;
        m->v[13] = y;
        m->v[14] = z;
}


void
mat4_vec4_mul(vec4_t *restrict res, mat4_t const *m, vec4_t const *restrict v)
{
        *res = VEC4(0.0, 0.0, 0.0, 0.0);

        for (int i = 0; i != 4; ++i) {
                float const *col = m->c[i];
                res->x += col[0]*v->v[i];
                res->y += col[1]*v->v[i];
                res->z += col[2]*v->v[i];
                res->w += col[3]*v->v[i];
        }
}


void
mat4_mul(mat4_t *res, mat4_t const *m0, mat4_t const *m1)
{
        mat4_t tmp;
        mat4_vec4_mul(&tmp.v0, m0, &m1->v0);
        mat4_vec4_mul(&tmp.v1, m0, &m1->v1);
        mat4_vec4_mul(&tmp.v2, m0, &m1->v2);
        mat4_vec4_mul(&tmp.v3, m0, &m1->v3);
        *res = tmp;
}


void
quat_mul(quat_t *res, quat_t const *restrict q0, quat_t const *restrict q1)
{
        res->w = q0->w*q1->w - q0->x*q1->x - q0->y*q1->y - q0->z*q1->z;
        res->x = q0->w*q1->x + q0->x*q1->w + q0->y*q1->z - q0->z*q1->y;
        res->y = q0->w*q1->y + q0->y*q1->w + q0->z*q1->x - q0->x*q1->z;
        res->z = q0->w*q1->z + q0->z*q1->w + q0->x*q1->y - q0->y*q1->x;
}


float
quat_dot(quat_t const *restrict q0, quat_t const *restrict q1)
{
        return q0->w*q1->w + q0->x*q1->x + q0->y*q1->y + q0->z*q1->z;
}


void
quat_normalize(quat_t *q)
{
        float const magnitude_inv = 1.0f/sqrtf(quat_dot(q, q));
        q->w *= magnitude_inv;
        q->x *= magnitude_inv;
        q->y *= magnitude_inv;
        q->z *= magnitude_inv;
}


void
quat_mat4(mat4_t *restrict m, quat_t const *restrict q)
{
        quat_t res = *q;
        /* check if quaternion is unit length */
        if (quat_dot(&res, &res) != 1.0f)
                quat_normalize(&res);

        m->v0 = VEC4(1.0f - 2.0f*(res.y*res.y + res.z*res.z),
                     2.0f*(res.x*res.y - res.w*res.z),
                     2.0f*(res.x*res.z + res.y*res.w), 0.0f);
        m->v1 = VEC4(2.0f*(res.x*res.y + res.w*res.z),
                     1.0f - 2.0f*(res.z*res.z + res.x*res.x),
                     2.0f*(res.z*res.y - res.w*res.x), 0.0f);
        m->v2 = VEC4(2.0f*(res.z*res.x - res.y*res.w),
                     2.0f*(res.z*res.y + res.w*res.x),
                     1.0f - 2.0f*(res.x*res.x + res.y*res.y), 0.0f);
        m->v3 = VEC4(0.0f, 0.0f, 0.0f, 1.0f);
}


void
quat_lerp(quat_t *restrict res,
          quat_t const *restrict q0,
          quat_t const *restrict q1,
          float const t)
{
        res->w = (1.0f - t)*q0->w + t*q1->w;
        res->x = (1.0f - t)*q0->x + t*q1->x;
        res->x = (1.0f - t)*q0->y + t*q1->y;
        res->z = (1.0f - t)*q0->z + t*q1->z;
        quat_normalize(res);
}


void
quat_slerp(quat_t *restrict res,
           quat_t const *restrict q0,
           quat_t const *restrict q1,
           float const t)
{
        float const dot = quat_dot(q0, q1);
        if (dot >= 1.0f) {
                *res = *q0;
                return;
        }

        /* small change therefore use linear interpolation */
        float const sin_omega = sqrtf(1.0f - dot*dot);
        if (sin_omega < 0.001f) {
                quat_lerp(res, q0, q1, t);
                return;
        }

        float const omega = acosf(dot);
        float const a = sinf((1.0f - t)*omega)/sin_omega;
        float const b = sinf(t*omega)/sin_omega;

        res->w = a*q0->w + b*q1->w;
        res->x = a*q0->x + b*q1->x;
        res->y = a*q0->y + b*q1->y;
        res->z = a*q0->z + b*q1->z;
}

#endif

#endif
