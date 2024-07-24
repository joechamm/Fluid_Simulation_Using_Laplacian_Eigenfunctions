#version 420

layout (location = 0) in int inIndex;

out int vIdx;

void main(void)
{
	vIdx = inIndex;
}