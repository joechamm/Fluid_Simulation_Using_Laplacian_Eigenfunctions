#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <glFramework/include/debug.h>

using glm::mat4;
using glm::vec4;
using glm::vec3;
using glm::vec2;

namespace glFramework {
	class GLFWApp
	{
		GLFWwindow* window_ = nullptr;
		double timeStamp_ = glfwGetTime();
		float deltaSeconds_ = 0;

	public:
		GLFWApp(int32_t width = 1600, int32_t height = 900, const char* title = "GLFW App Example")
		{
			glfwSetErrorCallback(
				[](int error, const char* description)
				{
					fprintf(stderr, "Error: %s\n", description);
				}
			);

			if (!glfwInit())
				exit(EXIT_FAILURE);

			glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
			glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
			glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifndef NDEBUG
			glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
#else
			glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_FALSE);
#endif

			if (width <= 0 || height <= 0)
			{
				const GLFWvidmode* info = glfwGetVideoMode(glfwGetPrimaryMonitor());
				window_ = glfwCreateWindow(info->width, info->height, (title ? title : "Full Screen Example"), nullptr, nullptr);
			}
			else
			{
				window_ = glfwCreateWindow(width, height, (title ? title : "Windowed Example"), nullptr, nullptr);
			}

			if (!window_)
			{
				glfwTerminate();
				exit(EXIT_FAILURE);
			}

			glfwMakeContextCurrent(window_);
			//gladLoadGL(glfwGetProcAddress);
			gladLoadGL();

			glfwSwapInterval(1);

#ifndef NDEBUG
			initDebug();
#endif

		}

		~GLFWApp()
		{
			glfwDestroyWindow(window_);
			glfwTerminate();
		}

		GLFWwindow* getWindow() const
		{
			return window_;
		}

		float getDeltaSeconds() const
		{
			return deltaSeconds_;
		}

		void swapBuffers()
		{
			glfwSwapBuffers(window_);
			glfwPollEvents();
			assert(glGetError() == GL_NO_ERROR);

			const double newTimeStamp = glfwGetTime();
			deltaSeconds_ = static_cast<float>(newTimeStamp - timeStamp_);
			timeStamp_ = newTimeStamp;
		}
	};
}