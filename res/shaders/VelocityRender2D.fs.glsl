#version 420

layout (location = 0) out vec4 FragColor;

uniform vec3 LineColor;

void main(void)
{
	FragColor = vec4(LineColor, 1.0);
}