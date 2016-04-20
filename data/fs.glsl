#version 450


out vec4 color;


in vec4 v_pos;


void
main(void)
{
        const vec3 light_dir = normalize(vec3(1, 0.5, 1));
        vec3 n = normalize(cross(dFdx(v_pos.xyz), dFdy(v_pos.xyz)));
        float geom_term = dot(n, light_dir);
        if (geom_term <= 0.0){
                color = vec4(0.2);
        } else {
                color = vec4(0.2) + vec4(0.7) * geom_term;
        }
}
