//#include "../include/GLShader.h"
//#include "../include/Utils.h"

#include <glFramework/include/GLShader.h>
#include <commonFramework/include/Utils.h>

#include <glad/glad.h>
#include <assert.h>
#include <stdio.h>
#include <string>

GLShader::GLShader(const char* filename)
	: GLShader(GLShaderTypeFromFilename(filename), readShaderFile(filename).c_str(), filename) {}

GLShader::GLShader(const std::string& filename)
	: GLShader(GLShaderTypeFromFilename(filename.c_str()), readShaderFile(filename.c_str()).c_str(), filename.c_str()) {}

GLShader::GLShader(GLenum type, const char* text, const char* debugFilename)
	: type_(type)
	, handle_(glCreateShader(type))
{
	glShaderSource(handle_, 1, &text, nullptr);
	glCompileShader(handle_);

	char buffer[8192];
	GLsizei length = 0;
	glGetShaderInfoLog(handle_, sizeof(buffer), &length, buffer);

	if (length)
	{
		printf("%s (File: %s)\n", buffer, debugFilename);
		printShaderSource(text);
		assert(false);
	}
}

GLShader::~GLShader()
{
	glDeleteShader(handle_);
}

void printProgramInfoLog(GLuint handle)
{
	char buffer[8192];
	GLsizei length = 0;
	glGetProgramInfoLog(handle, sizeof(buffer), &length, buffer);
	if (length)
	{
		printf("%s\n", buffer);
		assert(false);
	}
}

GLProgram::GLProgram(const GLShader& a)
	: handle_(glCreateProgram())
{
	glAttachShader(handle_, a.getHandle());
	glLinkProgram(handle_);
	printProgramInfoLog(handle_);
}

GLProgram::GLProgram(const GLShader& a, const GLShader& b)
	: handle_(glCreateProgram())
{
	glAttachShader(handle_, a.getHandle());
	glAttachShader(handle_, b.getHandle());
	glLinkProgram(handle_);
	printProgramInfoLog(handle_);
}

GLProgram::GLProgram(const GLShader& a, const GLShader& b, const GLShader& c)
	: handle_(glCreateProgram())
{
	glAttachShader(handle_, a.getHandle());
	glAttachShader(handle_, b.getHandle());
	glAttachShader(handle_, c.getHandle());
	glLinkProgram(handle_);
	printProgramInfoLog(handle_);
}

GLProgram::GLProgram(const GLShader& a, const GLShader& b, const GLShader& c, const GLShader& d, const GLShader& e)
	: handle_(glCreateProgram())
{
	glAttachShader(handle_, a.getHandle());
	glAttachShader(handle_, b.getHandle());
	glAttachShader(handle_, c.getHandle());
	glAttachShader(handle_, d.getHandle());
	glAttachShader(handle_, e.getHandle());
	glLinkProgram(handle_);
	printProgramInfoLog(handle_);
}

GLProgram::~GLProgram()
{
	glDeleteProgram(handle_);
}

void GLProgram::useProgram() const
{
	glUseProgram(handle_);
}

GLenum GLShaderTypeFromFileName(const char* filename)
{
	if (endsWith(filename, ".vert"))
		return GL_VERTEX_SHADER;

	if (endsWith(filename, ".frag"))
		return GL_FRAGMENT_SHADER;

	if (endsWith(filename, ".geom"))
		return GL_GEOMETRY_SHADER;

	if(endsWith(filename, ".tesc"))
		return GL_TESS_CONTROL_SHADER;

	if (endsWith(filename, ".tese"))
		return GL_TESS_EVALUATION_SHADER;

	if (endsWith(filename, ".comp"))
		return GL_COMPUTE_SHADER;

	assert(false);

	return 0;
}

GLBuffer::GLBuffer(GLsizeiptr size, const void* data, GLbitfield flags)
{
	glCreateBuffers(1, &handle_);
	glNamedBufferStorage(handle_, size, data, flags);
}

GLBuffer::~GLBuffer()
{
	glDeleteBuffers(1, &handle_);
}