#version 450


layout(location = 0) in vec3 pos;
layout(binding = 0, std140) uniform MatrixBlock {
        mat4 model;
        mat4 view;
        mat4 projection;
};


void
main(void)
{
        const mat4 model = mat4(
                 0.5,  0.0,  0.0,  0.0,
                 0.0,  0.0, -0.5,  0.0,
                 0.0,  0.5,  0.0,  0.0,
                 0.0,  0.0,  0.0,  1.0
        );

        gl_Position = projection*view*model*vec4(pos, 1.0);
}
