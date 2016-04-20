#version 450

layout(location = 0) uniform vec4 diffuse;

out vec4 color;


in vec4 v_pos;


void
main(void)
{
        vec3 n = normalize(cross(dFdx(v_pos.xyz), dFdy(v_pos.xyz)));
        color = diffuse*vec4(0.5*n + 0.5, 1.0);
//        color = vec4(1 - gl_FragCoord.z);
}
