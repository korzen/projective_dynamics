#version 330


layout(location = 0) in vec3 pos;


void
main(void)
{
        const mat4 model = mat4(
                 0.5,  0.0,  0.0,  0.0,
                 0.0,  0.0, -0.5,  0.0,
                 0.0,  0.5,  0.0,  0.0,
                 0.0,  0.0,  0.0,  1.0
        );

        gl_Position = model * vec4(pos, 1.0);
}
