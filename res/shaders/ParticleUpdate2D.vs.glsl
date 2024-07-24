#version 420

layout (location = 0) in vec2 CurrentPos2D;

layout (location = 0) out vec2 NewPos2D;

uniform samplerBuffer VelocityField;

uniform mat3 DomainToUnitTransform;
uniform mat3 UnitToDomainTransform;
uniform mat3 VelocityToUnitTransform;
uniform float dt;
uniform int GridWidth;
uniform int GridHeight;

void SampleVelocity(in vec2 normed_pos, out vec2 samp_vel)
{
	vec2 pos;
	vec2 vel[4];
	float weight[4];
	int idx[4];
	int col[2];
	int row[2];

	pos = normed_pos * vec2(float(GridWidth - 1), float(GridHeight - 1));

	col[0] = clamp(int(floor(pos.x)), 0, GridWidth);
	row[0] = clamp(int(floor(pos.y)), 0, GridHeight);
	col[1] = clamp((col[0] + 1), 0, GridWidth);
	row[1] = clamp((row[0] + 1), 0, GridHeight);

	idx[0] = (row[0] * GridWidth + col[0]);
	idx[1] = (row[1] * GridWidth + col[0]);
	idx[2] = (row[0] * GridWidth + col[1]);
	idx[3] = (row[1] * GridWidth + col[1]);

	/*
	 make the weight correspond to (1 - width) * (1 - height) (i.e. the larger the
	 area of the axis aligned rectangle defined by the corresponding grid point and pos,
	 the smaller the weight)
	*/

	float fracX = pos.x - float(col[0]);
	float fracY = pos.y - float(row[0]);

	if(col[0] == col[1])
	{
		fracX = 0.0;
	}
	if(row[0] == row[1])
	{
		fracY = 0.0;
	}
	
	weight[0] = (1.0 - fracX) * (1.0 - fracY);
	weight[1] = (1.0 - fracX) * fracY;
	weight[2] = fracX * (1.0 - fracY);
	weight[3] = fracX * fracY;
	
	vel[0] = texelFetch(VelocityField, idx[0]).xy;
	vel[1] = texelFetch(VelocityField, idx[1]).xy;
	vel[2] = texelFetch(VelocityField, idx[2]).xy;
	vel[3] = texelFetch(VelocityField, idx[3]).xy;

	vec2 avgVel = (vel[0] * weight[0] + vel[1] * weight[1] + vel[2] * weight[2] + vel[3] * weight[3]);

	samp_vel = vec2(VelocityToUnitTransform * vec3(avgVel, 0.0));
}

void main(void)
{
	vec2 unitPos, vel;
	unitPos = vec2(DomainToUnitTransform * vec3(CurrentPos2D, 1.0));
	SampleVelocity(unitPos, vel);
	unitPos = clamp(unitPos + vel * dt, 0.01, 0.99);
	NewPos2D = vec2(UnitToDomainTransform * vec3(unitPos, 1.0));
}