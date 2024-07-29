#pragma once

#include <glad/glad.h>
#include <string>

/// <summary>
///  The GLShader, GLProgram, and GLBuffer classes use the RAII idiom
/// </summary>
class GLShader
{
	GLenum type_;
	GLuint handle_;

public:
	explicit GLShader(const char* filename);
	GLShader(const std::string& filename);
	GLShader(GLenum type, const char* text, const char* debugFilename = "");
	~GLShader();

	GLenum getType() const
	{
		return type_;
	}

	GLuint getHandle() const
	{
		return handle_;
	}
};

class GLProgram
{
	GLuint handle_;
public:
	GLProgram(const GLShader& a);
	GLProgram(const GLShader& a, const GLShader& b);
	GLProgram(const GLShader& a, const GLShader& b, const GLShader& c);
	GLProgram(const GLShader& a, const GLShader& b, const GLShader& c, const GLShader& d, const GLShader& e);
	~GLProgram();

	void useProgram() const;

	GLuint getHandle() const
	{
		return handle_;
	}

};

/// return the shader type based on file's extension
GLenum GLShaderTypeFromFilename(const char* filename);

class GLBuffer
{
	GLuint handle_;
public:
	GLBuffer(GLsizeiptr size, const void* data, GLbitfield flags);
	~GLBuffer();

	GLuint getHandle() const
	{
		return handle_;
	}
};