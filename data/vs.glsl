#version 450


layout(location = 0) in vec3 pos;
layout(binding = 0, std140) uniform MatrixBlock {
        mat4 model;
        mat4 view;
        mat4 projection;
};


out vec4 v_pos;


void
main(void)
{
        v_pos = view*model*vec4(pos, 1.0);
        gl_Position = projection*v_pos;
}
