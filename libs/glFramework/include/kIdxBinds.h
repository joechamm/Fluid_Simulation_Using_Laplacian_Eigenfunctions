#pragma once

#include <glad/glad.h>
#include <cstdint>

// This file is used to define the binding points for the shader uniforms

namespace glFramework
{
	constexpr GLuint kIdxBind_PerFrameData = 7;
	constexpr GLuint kIdxBind_Materials = 6;
	constexpr GLuint kIdxBind_ModelMatrices = 5;
	constexpr GLuint kIdxBind_Vertices = 1;
}