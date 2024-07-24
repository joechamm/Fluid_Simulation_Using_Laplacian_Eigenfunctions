#version 420

layout (points) in;
layout (line_strip, max_vertices = 2) out;

in vec2 vPos[1];
in vec2 vVel[1];

uniform mat4 MVP;
uniform vec2 VelocityLengthScale;

void main(void)
{
	vec2 pt2 = vPos[0] + vVel[0] * VelocityLengthScale;

	gl_Position = MVP * vec4(vPos[0], 0.0, 1.0);
	EmitVertex();
	gl_Position = MVP * vec4(pt2, 0.0, 1.0);
	EmitVertex();
	EndPrimitive();
}