#pragma once

#include <commonFramework/include/scene/vec4.h>

namespace commonFramework
{

	// Common definitions for Lights and Cameras

	struct PACKED_STRUCT CameraProperties
	{
		gpumat4 projection_ = gpumat4(mat4(1.0f));
		gpumat4 view_ = gpumat4(mat4(1.0f));
		gpumat4 model_ = gpumat4(mat4(1.0f));
		gpuvec4 position_ = gpuvec4(vec4(1.0f));
	};

	static_assert(sizeof(CameraProperties) == (3 * 4 * sizeof(gpuvec4) + sizeof(gpuvec4)), "Invalid sizeof(CameraProperties), should be 3 * sizeof(mat4) + sizeof(gpuvec4)");

	struct PACKED_STRUCT LightProperties
	{
		gpuvec4 ambientColor_ = gpuvec4(0.0f);
		gpuvec4 color_ = gpuvec4(0.0f);
		// x = light radius, for point and spot lights; x = attenuation cutoff; 2 flag list floats [scattering + cast shadows]
		gpuvec4 attenuationParams_ = gpuvec4(1.0f, 2.0f / 255.0f, 0.0f, 0.0f);
		gpuvec4 FDirection = gpuvec4(0.0f, 0.0f, -1.0f, 0.0f);
		// inner angle, outer angle, falloff
		gpuvec4 FSpotRange = gpuvec4(0.0f, 90.0f, 1.0f, 0.0f);
		// (mapSize, mapCount, splitLambdaParam, strenght)
		gpuvec4 FShadowMapSize = gpuvec4(1024.0f, 0.0f, 0.5f, 1.0f);
		// const, slope, 0, 0 [ or filterType + shadowUpdateInterval]
		gpuvec4 shadowBias_ = gpuvec4(-0.0001f, 0.0f, 0.0f, 0.0f);
	};

	static_assert(sizeof(LightProperties) == 7 * sizeof(gpuvec4), "Invalid sizeof(LightProperties), should be 7 * sizeof(gpuvec4)");
}