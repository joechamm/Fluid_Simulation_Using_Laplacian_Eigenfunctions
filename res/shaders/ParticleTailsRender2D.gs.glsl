#version 420

layout (points) in;
layout (line_strip, max_vertices = 128) out;

in int vIdx[1];

out vec3 gColor;

layout (binding = 0) uniform samplerBuffer StoredPositions2DBuffer;

uniform mat4 MVP;
uniform vec3 HeadLineColor;
uniform vec3 TailLineColor;
uniform int NumParticles;
uniform int NumPositionsStored;
uniform int CurrentStoreOffset;

void main(void)
{
	int num_indices, idx;
	vec2 pos;
	float m, dm;
	dm = 1.0 / float(NumPositionsStored);
	num_indices = NumParticles * NumPositionsStored;
	idx = NumParticles * CurrentStoreOffset + vIdx[0];

	for(int i = 0; i < NumPositionsStored; i++)
	{
		m = float(i) * dm;
		gColor = mix(HeadLineColor, TailLineColor, m);
		pos = texelFetch(StoredPositions2DBuffer, idx).xy;
		gl_Position = MVP * vec4(pos, 0.0, 1.0);
		EmitVertex();

		idx = idx - NumParticles;
		if(idx < 0)
		{
			idx = num_indices + idx;
		}
	}
	
	EndPrimitive();
}