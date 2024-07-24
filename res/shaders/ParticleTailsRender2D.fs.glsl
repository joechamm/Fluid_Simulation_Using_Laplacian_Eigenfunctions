#version 420

layout (location = 0) out vec4 FragColor;

in vec3 gColor;

void main(void)
{
	FragColor = vec4(gColor, 1.0);
}