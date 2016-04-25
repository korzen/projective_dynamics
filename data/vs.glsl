#version 450


layout(location = 0) in vec3 pos;
layout(binding = 0, std140) uniform MatrixBlock {
        mat4 model;
        mat4 view;
        mat4 projection;
};
layout(location = 0) uniform int block_m;


out vec4 v_pos;
out vec4 v_color;

out gl_PerVertex {
        vec4 gl_Position;
        float gl_PointSize;
        float gl_ClipDistance[];
};

const vec4 colors[3] = vec4[](
        vec4(1.0, 0.0, 0.0, 1.0),
        vec4(0.0, 1.0, 0.0, 1.0),
        vec4(0.0, 0.0, 1.0, 0.0)
);


void
main(void)
{
        /* matrix row contains all positions, so we can show only row blocks */
        int block_x = gl_VertexID/block_m;
        v_color = colors[block_x%3];

        v_pos = view*model*vec4(pos, 1.0);
        gl_Position = projection*v_pos;
}
