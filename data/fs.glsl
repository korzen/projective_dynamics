#version 450


out vec4 color;


in vec4 v_pos;


void
main(void)
{
        vec3 n = normalize(cross(dFdx(v_pos.xyz), dFdy(v_pos.xyz)));
        color = vec4(0.5*n + 0.5, 1.0);
//        color = vec4(1 - gl_FragCoord.z);
}
