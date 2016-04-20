#version 450

layout(location = 0) uniform vec4 diffuse;

out vec4 color;


in vec4 v_pos;


void
main(void)
{
        color = vec4(diffuse);
}
