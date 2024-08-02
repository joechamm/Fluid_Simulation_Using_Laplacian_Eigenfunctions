#include <helpers/RootDir.h>
#include <glFramework/include/GLFWApp.h>
#include <glFramework/include/GLShader.h>

#include <commonFramework/include/ResourceString.h>
#include <commonFramework/include/UtilsFPS.h>

#include <fluidSim/include/FluidSimulation.h>
#include <fluidSim/include/VelocityField.h>

#include <imgui.h>

using glm::mat4;
using glm::vec4;
using glm::dvec4;
using glm::dvec2;
using glm::vec2;
using glm::ivec2;

struct PerFrameData
{
	mat4 view;
	mat4 proj;
	vec4 cameraPos;
};

struct MouseState
{
	vec2 pos = vec2(0.0f);
	bool pressedLeft = false;
} mouseState;

GLFWwindow* g_window = nullptr;

const ivec2 g_gridSize = ivec2(36, 36);
const double g_viscosity = 0.01;
const double g_dt = 0.1;
const double g_amplitudeScale = 1.0;

const float g_velocityLengthScale = 30.0f;

const vec3 g_lineColor = vec3(0.0f, 0.8f, 0.4f);

GLuint g_vao = 0;

void InitGridBuffer(GLuint gridBufferHandle, int gridSizeX, int gridSizeY);
void ExpandBasis(GLuint velocityBufferHandle, FluidSimulation& sim, VelocityField& field);


int main(int argc, char** argv)
{
	glFramework::GLFWApp app(1600, 900, "Fluid Simulation");

	g_window = app.getWindow();

	PerFrameData perFrameData;
	perFrameData.proj = glm::ortho(0.0f, 1600.0f, 0.0f, 900.0f);

	commonFramework::FramesPerSecondCounter fpsCounter;

	FluidSimulation fluidSim(6, 0.01, 0.1, dvec4(0.0, 0.0, M_PI, M_PI));
	VelocityField velocityField(g_gridSize.x, g_gridSize.y, 6, g_amplitudeScale);

	glFramework::GLBuffer gridBuffer(g_gridSize.x * g_gridSize.y * sizeof(GLfloat), nullptr, GL_DYNAMIC_STORAGE_BIT);
	glFramework::GLBuffer velocityBuffer(g_gridSize.x * g_gridSize.y * sizeof(GLfloat), nullptr, GL_DYNAMIC_STORAGE_BIT);
	
	glCreateVertexArrays(1, &g_vao);
	glBindVertexArray(g_vao);
	glVertexArrayVertexBuffer(g_vao, 0, gridBuffer.getHandle(), 0, 2 * sizeof(GLfloat));
	glVertexArrayVertexBuffer(g_vao, 1, velocityBuffer.getHandle(), 0, 2 * sizeof(GLfloat));

	fluidSim.Initialize();
	velocityField.PrecomputeBasisField(fluidSim);

	// Set up the shader
	glFramework::GLShader vtxVelocity(commonFramework::appendToRoot("res/shaders/velocity.vs"));
	glFramework::GLShader geomVelocity(commonFramework::appendToRoot("res/shaders/velocity.gs"));
	glFramework::GLShader fragVelocity(commonFramework::appendToRoot("res/shaders/velocity.fs"));
	glFramework::GLProgram velocityProgram(vtxVelocity, fragVelocity, geomVelocity);
	velocityProgram.useProgram();

	glfwSetCursorPosCallback(
		app.getWindow(),
		[](auto* window, double x, double y)
		{
			int width, height;
			glfwGetFramebufferSize(window, &width, &height);
			mouseState.pos.x = static_cast<float>(x / width);
			mouseState.pos.y = static_cast<float>(y / height);
			ImGui::GetIO().MousePos = ImVec2((float)x, (float)y);
		}
	);

	glfwSetMouseButtonCallback(
		app.getWindow(),
		[](auto* window, int button, int action, int mods)
		{
			auto& io = ImGui::GetIO();
			const int idx = button == GLFW_MOUSE_BUTTON_LEFT ? 0 : button == GLFW_MOUSE_BUTTON_RIGHT ? 2 : 1;
			io.MouseDown[idx] = action == GLFW_PRESS;

			if (!io.WantCaptureMouse)
				if (button == GLFW_MOUSE_BUTTON_LEFT)
					mouseState.pressedLeft = action == GLFW_PRESS;
		}
	);

	glfwSetKeyCallback(
		app.getWindow(),
		[](GLFWwindow* window, int key, int scancode, int action, int mods)
		{
			const bool pressed = action != GLFW_RELEASE;
			if (key == GLFW_KEY_ESCAPE && pressed)
				glfwSetWindowShouldClose(window, GLFW_TRUE);

		}
	);

	// Set up the mouse callback
	glfwSetCursorPosCallback(app.getWindow(), [](GLFWwindow* window, double xpos, double ypos)
		{
			mouseState.pos = vec2(xpos, ypos);
		});

	glfwSetMouseButtonCallback(app.getWindow(), [](GLFWwindow* window, int button, int action, int mods)
		{
			if (button == GLFW_MOUSE_BUTTON_LEFT)
			{
				if (action == GLFW_PRESS)
				{
					mouseState.pressedLeft = true;
				}
				else if (action == GLFW_RELEASE)
				{
					mouseState.pressedLeft = false;
				}
			}
		});

	// Main loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{
		/// Update
		fpsCounter.tick(app.getDeltaSeconds());
		fluidSim.TimeStep();
		ExpandBasis(velocityBuffer.getHandle(), fluidSim, velocityField);
		
		/// Render
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		velocityProgram.useProgram();


		glUniformMatrix4fv(0, 1, GL_FALSE, glm::value_ptr(perFrameData.proj * perFrameData.view));
		glUniform2f(1, g_velocityLengthScale, g_velocityLengthScale);
		glUniform3f(2, g_lineColor.r, g_lineColor.g, g_lineColor.b);

		glBindVertexArray(g_vao);
		glDrawArrays(GL_POINTS, 0, g_gridSize.x * g_gridSize.y);

		app.swapBuffers();
	}

	return 0;
}



void InitGridBuffer(GLuint gridBufferHandle, int gridSizeX, int gridSizeY)
{
	std::vector<GLfloat> gridData;
	gridData.reserve(gridSizeX * gridSizeY * 2);

	for (int y = 0; y != gridSizeY; y++)
	{
		for (int x = 0; x != gridSizeX; x++)
		{
			gridData.push_back(static_cast<GLfloat>(x) / static_cast<GLfloat>(gridSizeX));
			gridData.push_back(static_cast<GLfloat>(y) / static_cast<GLfloat>(gridSizeY));
		}
	}

	glNamedBufferSubData(gridBufferHandle, 0, gridData.size() * sizeof(GLfloat), gridData.data());
}

void ExpandBasis(GLuint velocityBufferHandle, FluidSimulation& sim, VelocityField& field)
{
	std::vector<GLfloat> velocityData;
	velocityData.reserve(g_gridSize.x * g_gridSize.y * 2);

	for (int row = 0; row <= g_gridSize.y; row++)
	{
		for (int col = 0; col <= g_gridSize.x; col++)
		{
			dvec2 v = field.ComputeVelocity(row, col, sim);
			velocityData.push_back(static_cast<GLfloat>(v.x));
			velocityData.push_back(static_cast<GLfloat>(v.y));
		}
	}

	glNamedBufferSubData(velocityBufferHandle, 0, g_gridSize.x * g_gridSize.y * 2 * sizeof(GLfloat), velocityData.data());
}