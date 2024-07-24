#version 420

layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aVel;

out vec2 vPos;
out vec2 vVel;

void main(void)
{
	vPos = aPos;
	vVel = aVel;
}