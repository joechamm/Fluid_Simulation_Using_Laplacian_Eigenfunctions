#include <iostream>
#include <fstream>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <assert.h>
#include <time.h>
#include <vector>
#include <armadillo>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/noise.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace std;
using glm::mat4;
using glm::mat3;
using glm::mat2;
using glm::vec4;
using glm::vec3;
using glm::vec2;

using glm::dmat4;
using glm::dmat3;
using glm::dmat2;
using glm::dvec4;
using glm::dvec3;
using glm::dvec2;

using glm::ivec4;
using glm::ivec3;
using glm::ivec2;

#define WINDOW_TITLE_PREFIX "Eigen Fluid Sim"

static vec3 VelocityColor(0.0f, 1.0f, 0.0f);

static mat4 ParticlesMVP;
static mat4 RenderMVP;

static mat3 UpdateTransform;
static mat3 UpdateInverseTransform;

const double PI = 3.14159265358979323846264;
const double DOMAIN_MIN_X = 0.0;
const double DOMAIN_MIN_Y = 0.0;
const double DOMAIN_MAX_X = PI;
const double DOMAIN_MAX_Y = PI;
const double DOMAIN_WIDTH = PI;
const double DOMAIN_HEIGHT = PI;
const int DIM = 2;

static vec2 velocity_scales;
static float DOMAIN_MARGIN = 0.00001f;

static int N = 36;
static int N_SQ = 6;
static int N_OLD = N;
static int VEL_GRID_RES_X = 36;
static int VEL_GRID_RES_Y = 36;

const int NUM_PARTICLES = 1000;
const int NUM_POSITIONS_STORED = 50;
const GLsizeiptr BUFFER_UPDATE_SIZE = NUM_PARTICLES * DIM * sizeof(float);
const GLsizeiptr BUFFER_TOTAL_SIZE = BUFFER_UPDATE_SIZE * NUM_POSITIONS_STORED;
static GLsizeiptr CURRENT_STORE_OFFSET = 0;
static int CURRENT_INDEX_OFFSET = 0;

static GLuint CURRENT_MODE = 0;
static int* zero_modes = 0;
bool ZERO_MODES_ON = false;

static ivec2* basis_lookup_table = 0;
static int** basis_rlookup_table = 0;
static dvec2*** velocity_basis = 0;
static double* eigs = 0;
static double* eigs_inv = 0;
static double* eigs_inv_root = 0;

static arma::sp_mat* CkMatrices = NULL;
static arma::vec					dwCoefficients;
static arma::vec					wCoefficients;
static arma::vec					dwForceCoefficients;
static arma::vec					rk4Qn[4];
static arma::vec					rk4Dwt[4];

static double		viscosity = 0.01;
static double		dt = 0.1;
static double		initial_amp_scale = 1.0;

unsigned FrameCount = 0;

static int CurrentWidth = 800;
static int CurrentHeight = 600;
static int MainWindow = 0;

static bool UPDATE_VELOCITY_GRID_ON = false;
static bool UPDATE_PARTICLES_ON = false;

/*
enum RENDER_MODE
{
	RENDER_VELOCITY_GRID,
	RENDER_PARTICLES,
	RENDER_GRID_AND_PARTICLES,
	RENDER_NUM_MODES
} CURRENT_RENDER_MODE = RENDER_PARTICLES;
*/

enum RENDER_MODE
{
	RENDER_VELOCITY_GRID,
	RENDER_PARTICLES,
	RENDER_GRID_AND_PARTICLES,
	RENDER_NUM_MODES
};

static unsigned int CURRENT_RENDER_MODE = RENDER_PARTICLES;

static float vel_render_len_scale = 30.0f;

static std::vector<ivec2>	drag_path;

static bool left_mouse_down = false;

static vec3 particle_head_color(0.0f, 0.0f, 1.0f);
static vec3 particle_tail_color(1.0f, 0.0f, 0.0f);
static float particle_point_scale = 3.5f;
static float particle_tail_line_width = 3.5f;

// params
static int visc_coarse = 0;
static int visc_fine = 100;

// force mode
static int force_mode = 2;

// mouse force multiplier
static double force_mult = 250.0;

// force vector on mouse drag
static dvec2 mouse_p1(0.0, 0.0);
static dvec2 mouse_p2(0.0, 0.0);

// update stuff

static bool vel_buff_size_changed = false;
static bool grid_pos_needs_resize = true;
static bool grid_vel_needs_resize = true;

GLuint
particleUpdatesProgHandle,
particleTailsProgHandle,
velocitiesProgHandle;

GLuint
grid_vao,
grid_pos_vbo,
grid_vel_vbo;

GLuint
particles_render_vao,
particles_update_vao;

GLuint
particles_indices_vbo,
particles_position_vbo;

GLuint
particles_dest_vbo;

GLuint
particles_position_tex,
velocity_tex;

GLint
vel_mvp_loc,
vel_col_loc,
vel_len_scale_loc,
particles_mvp_loc,
particles_head_col_loc,
particles_tail_col_loc,
particles_num_particles_loc,
particles_num_stored_loc,
particles_current_store_loc,
particles_tex_buff_loc,
update_trans_loc,
update_inv_trans_loc,
update_vel_scale_loc,
update_dt_loc,
update_tex_loc,
update_width_loc,
update_height_loc;

// initialize fields and tables
void InitFields(void);

void CreateLookupTables(void);
void CreateCkMatrices(void);
void CreateBasisField(void);
void CreateEigenValues(void);

void FillLookupTable(void);
void PrecomputeBasisFields(void);
void PrecomputeDynamics(void);
void BasisFieldRect2D(const ivec2& K, double amp, dvec2** vBasisElem);

// destroy fields and tables
void DestroyLookupTables(void);
void DestroyCkMatrices(void);
void DestroyBasisField(void);
void DestroyEigenValues(void);

// energy functions
double		CurrentEnergy(void);
void			SetEnergy(double energy);

// coefficients
double CoefDensity(const ivec2& a, const ivec2& b, const ivec2& c);

// index lookup functions
int BasisLookup(int index, int component);
ivec2 BasisLookupK(int index);
int BasisRLookup(const ivec2& K);

// write velocity field as sum of basis fields
void ExpandBasis(void);

// force functions
void ProjectForces(const std::vector<dvec4>& force_path, arma::vec& forceVec);
void Stir(const std::vector<dvec4>& force_path);

// updates
void UpdateFields(void);

// initializations
void InitGlobals(void);
void InitBuffers(void);
void InitTextures(void);
void InitShaders(void);
void InitUniforms(void);
void Init(int, char**);

// loading shader from file
GLuint LoadShaderFromFile(GLenum, const char*);

// render sub functions
void RenderGrid(void);

// glut callbacks
void Render(void);
void Idle(void);
void Timer(int);
void Keyboard(unsigned char, int, int);
void SpecialKey(int, int, int);
void Reshape(int, int);
void Mouse(int, int, int, int);
void MouseMove(int, int);

// reset
void Reset(void);

// shutdown
void Shutdown(void);

// playing with the basis dimension
void SetBasisDimension(int dim);

// playing with velocity dimension
void SetVelocityResolution(int res);

void CreateGridRenderArray(void);
void DestroyGridRenderArray(void);

void InitGrid(void);
void DestroyGrid(void);

void CreateGridVelocityBuffer(void);
void CreateGridPositionBuffer(void);
void FillGridVelocityBuffer(void);
void FillGridPositionBuffer(void);
void DestroyGridVelocityBuffer(void);
void DestroyGridPositionBuffer(void);

void CreateParticleBuffers(void);
void FillParticleBuffers(void);
void DestroyParticleBuffers(void);

void CreateParticleRenderArrays(void);
void DestroyParticleRenderArrays(void);

void CreateParticleUpdateArrays(void);
void DestroyParticleUpdateArrays(void);

void RenderParticles(void);
void UpdateParticles(void);

void CopyParticleBuffers(void);

void CreateVelocityTexture(void);
void DestroyVelocityTexture(void);
void CreatePositionTexture(void);
void DestroyPositionTexture(void);

void SetRenderMode(int mode);

void ZeroModes(void);
void ZeroMultiplesOf(int m);
void InsertMultiplesOf(int m);

int main(int argc, char** argv)
{
	viscosity = (double(visc_coarse) + double(visc_fine) / 200.0) / 500.0;

	N_SQ = int(glm::floor(sqrt(double(N))));
	InitGlobals();
	srand(time(0));
	printf("Initializing Fields\n");
	InitFields();
	printf("Fill Lookups\n");
	FillLookupTable();
	printf("Precompute Basis Fields\n");
	PrecomputeBasisFields();
	printf("Precompute Dynamics\n");
	PrecomputeDynamics();
	printf("Init Uniforms\n");
	InitUniforms();
	printf("Initializing Glut\n");
	Init(argc, argv);
	printf("Initializing Shaders\n");
	InitShaders();
	printf("Initializing Buffers\n");
	InitBuffers();
	printf("Initializing Textures\n");
	InitTextures();

	glutMainLoop();
	Shutdown();
	return 0;
}

void SetRenderMode(int mode)
{
	if (mode < 0 || mode >= int(RENDER_NUM_MODES))
	{
		printf("Error: Mode Out of Range\n");
		return;
	}

	CURRENT_RENDER_MODE = (RENDER_MODE)mode;
}

void Render(void)
{
	FrameCount++;

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	switch (CURRENT_RENDER_MODE)
	{
	case RENDER_VELOCITY_GRID:
		RenderGrid();
		break;
	case RENDER_PARTICLES:
		RenderParticles();
		break;
	case RENDER_GRID_AND_PARTICLES:
		RenderGrid();
		RenderParticles();
		break;
	default:
		RenderParticles();
		RenderGrid();
		break;
	}

	glutSwapBuffers();
	glutPostRedisplay();
}

void RenderGrid(void)
{
	vec2 velocity_length_scale = vec2(vel_render_len_scale / float(DOMAIN_WIDTH), vel_render_len_scale / float(DOMAIN_HEIGHT));

	GLuint numGridPoints = (VEL_GRID_RES_X + 1) * (VEL_GRID_RES_Y + 1);

	glLineWidth(1.0f);

	glUseProgram(velocitiesProgHandle);

	glUniformMatrix4fv(vel_mvp_loc, 1, GL_FALSE, glm::value_ptr(RenderMVP));
	glUniform3fv(vel_col_loc, 1, glm::value_ptr(VelocityColor));
	glUniform2fv(vel_len_scale_loc, 1, glm::value_ptr(velocity_length_scale));

	glBindVertexArray(grid_vao);

	glDrawArrays(GL_POINTS, 0, numGridPoints);

	glBindVertexArray(0);

	glUseProgram(0);
}

void RenderParticles(void)
{
	ParticlesMVP = glm::ortho(float(DOMAIN_MIN_X), float(DOMAIN_MAX_X), float(DOMAIN_MIN_Y), float(DOMAIN_MAX_Y));

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_BUFFER, particles_position_tex);

	glLineWidth(particle_tail_line_width);

	glUseProgram(particleTailsProgHandle);

	glUniformMatrix4fv(particles_mvp_loc, 1, GL_FALSE, glm::value_ptr(ParticlesMVP));
	glUniform3fv(particles_head_col_loc, 1, glm::value_ptr(particle_head_color));
	glUniform3fv(particles_tail_col_loc, 1, glm::value_ptr(particle_tail_color));
	glUniform1i(particles_num_particles_loc, NUM_PARTICLES);
	glUniform1i(particles_num_stored_loc, NUM_POSITIONS_STORED);
	glUniform1i(particles_current_store_loc, CURRENT_INDEX_OFFSET);
	glUniform1i(particles_tex_buff_loc, 0);

	glBindVertexArray(particles_render_vao);

	glDrawArrays(GL_POINTS, 0, NUM_PARTICLES);

	glBindVertexArray(0);

	glUseProgram(0);

	glBindTexture(GL_TEXTURE_BUFFER, 0);
}

void Reset(void)
{
	UPDATE_VELOCITY_GRID_ON = false;
	UPDATE_PARTICLES_ON = false;

	printf("Resetting\n");
	InitFields();
	printf("Filling Lookup Tables\n");
	FillLookupTable();
	printf("Precomputing Basis Fields\n");
	PrecomputeBasisFields();
	printf("Precomputing Dynamics\n");
	PrecomputeDynamics();

	FillGridPositionBuffer();
	FillGridVelocityBuffer();
	DestroyVelocityTexture();
	CreateVelocityTexture();
	FillParticleBuffers();

	ExpandBasis();
}

void ZeroModes(void)
{
	for (int k = 0; k < N; k++)
	{
		if (zero_modes[k] == 0)
		{
			wCoefficients(k) = 0.0;
		}
	}
}

void ZeroMultiplesOf(int m)
{
	if (m == 0)
	{
		zero_modes[0] = 0;
		return;
	}

	for (int i = 1; i < N; i++)
	{
		if (i % m == 0)
		{
			zero_modes[i] = 0;
		}
	}
}

void InsertMultiplesOf(int m)
{
	if (m == 0)
	{
		zero_modes[0] = 1;
		return;
	}

	for (int i = 1; i < N; i++)
	{
		if (i % m == 0)
		{
			zero_modes[i] = 1;
		}
	}
}

void UpdateFields(void)
{
	double prev_e = CurrentEnergy();

	// 4th order Runge-Kutta to update velocity
	rk4Qn[0] = wCoefficients;

	for (int k = 0; k < N; k++)
	{
		rk4Dwt[0](k) = arma::as_scalar(arma::trans(rk4Qn[0]) * CkMatrices[k] * rk4Qn[0]);
		rk4Qn[1](k) = rk4Qn[0](k) + rk4Dwt[0](k) * dt * 0.5;
	}

	for (int k = 0; k < N; k++)
	{
		rk4Dwt[1](k) = arma::as_scalar(arma::trans(rk4Qn[1]) * CkMatrices[k] * rk4Qn[1]);
		rk4Qn[2](k) = rk4Qn[0](k) + rk4Dwt[1](k) * dt * 0.5;
	}

	for (int k = 0; k < N; k++)
	{
		rk4Dwt[2](k) = arma::as_scalar(arma::trans(rk4Qn[2]) * CkMatrices[k] * rk4Qn[2]);
		rk4Qn[3](k) = rk4Qn[0](k) + rk4Dwt[2](k) * dt * 0.5;
	}

	for (int k = 0; k < N; k++)
	{
		rk4Dwt[3](k) = arma::as_scalar(arma::trans(rk4Qn[3]) * CkMatrices[k] * rk4Qn[3]);
		dwCoefficients(k) = (rk4Dwt[0](k) + rk4Dwt[1](k) * 2.0 + rk4Dwt[2](k) * 2.0 + rk4Dwt[3](k)) / 6.0;
	}

	wCoefficients += dwCoefficients * dt;

	if (prev_e > 1e-5)
	{
		SetEnergy(prev_e);
	}


	for (int k = 0; k < N; k++)
	{
		double eig = -eigs[k];
		wCoefficients(k) = wCoefficients(k) * exp(eig * dt * viscosity) + dwForceCoefficients(k);
		dwForceCoefficients(k) = 0.0;
	}

	ExpandBasis();
}

void InitGlobals(void)
{
	velocitiesProgHandle = 0;
	particleTailsProgHandle = 0;
	particleUpdatesProgHandle = 0;

	grid_vao = 0;
	grid_pos_vbo = 0;
	grid_vel_vbo = 0;

	particles_render_vao = 0;
	particles_update_vao = 0;

	particles_indices_vbo = 0;
	particles_position_vbo = 0;
	particles_dest_vbo = 0;

	particles_position_tex = 0;
	velocity_tex = 0;

	vel_mvp_loc = -1;
	vel_col_loc = -1;
	vel_len_scale_loc = -1;
	particles_mvp_loc = -1;
	particles_head_col_loc = -1;
	particles_tail_col_loc = -1;
	particles_num_particles_loc = -1;
	particles_num_stored_loc = -1;
	particles_current_store_loc = -1;
	particles_tex_buff_loc = -1;
	update_trans_loc = -1;
	update_inv_trans_loc = -1;
	update_vel_scale_loc = -1;
	update_dt_loc = -1;
	update_tex_loc = -1;
	update_width_loc = -1;
	update_height_loc = -1;

	if (zero_modes != NULL)
	{
		delete[] zero_modes;
		zero_modes = NULL;
	}

	zero_modes = new int[N];

	for (int k = 0; k < N; k++)
	{
		zero_modes[k] = 1;
	}

	for (int i = 0; i < 4; i++)
	{
		rk4Qn[i].set_size(N);
		rk4Dwt[i].set_size(N);
	}
}

void InitBuffers(void)
{
	printf("Initializing Grid\n");
	InitGrid();
	printf("Creating Partricles Buffers\n");
	CreateParticleBuffers();
	printf("Filling Particles Buffers\n");
	FillParticleBuffers();
	printf("Createting Particles Render Arrays\n");
	CreateParticleRenderArrays();
	printf("Creating Updates\n");
	CreateParticleUpdateArrays();
}

void InitTextures(void)
{
	CreateVelocityTexture();
	CreatePositionTexture();
}

void CreateGridRenderArray(void)
{
	if (!glIsBuffer(grid_pos_vbo))
	{
		printf("grid_pos_vbo is NOT buffer object\n");
		CreateGridPositionBuffer();
		FillGridPositionBuffer();
	}
	if (!glIsBuffer(grid_vel_vbo))
	{
		printf("grid_vel_vbo is NOT buffer object\n");
		CreateGridVelocityBuffer();
		FillGridVelocityBuffer();
	}
	if (glIsVertexArray(grid_vao))
	{
		DestroyGridRenderArray();
	}

	glGenVertexArrays(1, &grid_vao);

	glBindVertexArray(grid_vao);

	glBindBuffer(GL_ARRAY_BUFFER, grid_pos_vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, grid_vel_vbo);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);

	glBindVertexArray(0);
}

void DestroyGridRenderArray(void)
{
	glBindVertexArray(0);
	if (glIsVertexArray(grid_vao))
	{
		glDeleteVertexArrays(1, &grid_vao);
		grid_vao = 0;
	}
}

void InitGrid(void)
{
	printf("Creating Grid Position Buffer\n");
	CreateGridPositionBuffer();
	printf("Filling Grid Position Buffer\n");
	FillGridPositionBuffer();
	printf("Creating Grid Velocity Buffer\n");
	CreateGridVelocityBuffer();
	printf("Filling Grid Velocity Buffer\n");
	FillGridVelocityBuffer();
	printf("Creating Grid Render Array\n");
	CreateGridRenderArray();
}

void DestroyGrid(void)
{
	DestroyGridVelocityBuffer();
	DestroyGridPositionBuffer();
	DestroyGridRenderArray();
}

void CreateGridVelocityBuffer(void)
{
	if (!glIsBuffer(grid_vel_vbo))
	{
		if (grid_vel_vbo)
		{
			glDeleteBuffers(1, &grid_vel_vbo);
			grid_vel_vbo = 0;
		}

		glGenBuffers(1, &grid_vel_vbo);
	}

	GLsizeiptr velBuffsize = (VEL_GRID_RES_X + 1) * (VEL_GRID_RES_Y + 1) * DIM * sizeof(float);

	glBindBuffer(GL_ARRAY_BUFFER, grid_vel_vbo);
	glBufferData(GL_ARRAY_BUFFER, velBuffsize, NULL, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	grid_vel_needs_resize = false;
}

void CreateGridPositionBuffer(void)
{
	if (!glIsBuffer(grid_pos_vbo))
	{
		if (grid_pos_vbo)
		{
			glDeleteBuffers(1, &grid_pos_vbo);
			grid_pos_vbo = 0;
		}

		glGenBuffers(1, &grid_pos_vbo);
	}

	GLsizeiptr posBuffsize = (VEL_GRID_RES_X + 1) * (VEL_GRID_RES_Y + 1) * DIM * sizeof(float);

	glBindBuffer(GL_ARRAY_BUFFER, grid_pos_vbo);
	glBufferData(GL_ARRAY_BUFFER, posBuffsize, NULL, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	grid_pos_needs_resize = false;
}

void DestroyGridVelocityBuffer(void)
{
	if (grid_vel_vbo)
	{
		glDeleteBuffers(1, &grid_vel_vbo);
		grid_vel_vbo = 0;
	}
}

void DestroyGridPositionBuffer(void)
{
	if (grid_pos_vbo)
	{
		glDeleteBuffers(1, &grid_pos_vbo);
		grid_pos_vbo = 0;
	}
}

void FillGridPositionBuffer(void)
{
	if (!glIsBuffer(grid_pos_vbo))
	{
		printf("grid_pos_vbo is NOT buffer object\n");
		CreateGridPositionBuffer();
	}

	float gridWidth = float(CurrentWidth);
	float gridHeight = float(CurrentHeight);
	float gridMinX = 0.0f;
	float gridMinY = 0.0f;
	float dx = gridWidth / float(VEL_GRID_RES_X);
	float dy = gridHeight / float(VEL_GRID_RES_Y);

	float* grid_positions = new float[(DIM * (VEL_GRID_RES_X + 1) * (VEL_GRID_RES_Y + 1))];

	int idx = 0;

	for (int row = 0; row <= VEL_GRID_RES_Y; row++)
	{
		float y = float(row) * dy + gridMinY;
		for (int col = 0; col <= VEL_GRID_RES_X; col++)
		{
			float x = float(col) * dx + gridMinX;
			grid_positions[idx] = x;
			idx++;
			grid_positions[idx] = y;
			idx++;
		}
	}

	GLsizeiptr gridPosSize = ((VEL_GRID_RES_X + 1) * (VEL_GRID_RES_Y + 1) * DIM * sizeof(float));

	glBindBuffer(GL_ARRAY_BUFFER, grid_pos_vbo);

	if (grid_pos_needs_resize)
	{
		glBufferData(GL_ARRAY_BUFFER, gridPosSize, NULL, GL_STATIC_DRAW);
		grid_pos_needs_resize = false;
	}

	glBufferSubData(GL_ARRAY_BUFFER, 0, gridPosSize, grid_positions);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	delete[] grid_positions;
}

void FillGridVelocityBuffer(void)
{
	if (!glIsBuffer(grid_vel_vbo))
	{
		printf("grid_vel_vbo is NOT buffer object\n");
		CreateGridVelocityBuffer();
	}

	GLsizeiptr gridVelSize = ((VEL_GRID_RES_X + 1) * (VEL_GRID_RES_Y + 1) * DIM * sizeof(float));

	printf("creating tmp vel array\n");
	float* vel_field = (float*)malloc(gridVelSize);
	printf("setting tmp vel array to 0\n");
	memset(vel_field, 0, gridVelSize);
	printf("binding buffer\n");
	glBindBuffer(GL_ARRAY_BUFFER, grid_vel_vbo);

	if (grid_vel_needs_resize)
	{
		printf("setting buffer data\n");
		glBufferData(GL_ARRAY_BUFFER, gridVelSize, NULL, GL_DYNAMIC_DRAW);
		grid_vel_needs_resize = false;
	}

	printf("writing buffer data\n");
	glBufferSubData(GL_ARRAY_BUFFER, 0, gridVelSize, vel_field);

	printf("unbinding vel buffer\n");
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	printf("deleting tmp vel array\n");
	free(vel_field);
	vel_field = 0;
}

void CreateParticleRenderArrays(void)
{
	glBindVertexArray(0);

	if (particles_render_vao != 0)
	{
		glDeleteVertexArrays(1, &particles_render_vao);
		particles_render_vao = 0;
	}

	if (!glIsBuffer(particles_indices_vbo) || !glIsBuffer(particles_position_vbo))
	{
		CreateParticleBuffers();
	}

	glGenVertexArrays(1, &particles_render_vao);

	glBindVertexArray(particles_render_vao);

	glBindBuffer(GL_ARRAY_BUFFER, particles_indices_vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribIPointer(0, 1, GL_INT, 0, 0);

	glBindVertexArray(0);
}

void DestroyParticleRenderArrays(void)
{
	glBindVertexArray(0);
	if (particles_render_vao != 0)
	{
		glDeleteVertexArrays(1, &particles_render_vao);
		particles_render_vao = 0;
	}
}

void CreateParticleBuffers(void)
{
	glBindVertexArray(0);

	if (particles_position_vbo != 0)
	{
		glDeleteBuffers(1, &particles_position_vbo);
		particles_position_vbo = 0;
	}

	if (particles_dest_vbo != NULL)
	{
		glDeleteBuffers(1, &particles_dest_vbo);
		particles_dest_vbo = 0;
	}

	if (particles_indices_vbo != NULL)
	{
		glDeleteBuffers(1, &particles_indices_vbo);
		particles_indices_vbo = 0;
	}

	GLsizeiptr indBuffSize = NUM_PARTICLES * sizeof(int);

	glGenBuffers(1, &particles_position_vbo);
	glGenBuffers(1, &particles_indices_vbo);
	glGenBuffers(1, &particles_dest_vbo);

	glBindBuffer(GL_ARRAY_BUFFER, particles_indices_vbo);
	glBufferData(GL_ARRAY_BUFFER, indBuffSize, NULL, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_ARRAY_BUFFER, particles_position_vbo);
	glBufferData(GL_ARRAY_BUFFER, BUFFER_TOTAL_SIZE, NULL, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, particles_dest_vbo);
	glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, BUFFER_UPDATE_SIZE, NULL, GL_DYNAMIC_COPY);
	glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);
}

void FillParticleBuffers(void)
{
	GLsizeiptr buff_size = NUM_PARTICLES * sizeof(int);

	if (!glIsBuffer(particles_position_vbo) || !glIsBuffer(particles_indices_vbo) || !glIsBuffer(particles_dest_vbo))
	{
		CreateParticleBuffers();
	}

	vec2* particle_positions = new vec2[NUM_PARTICLES];
	int* particle_indices = new int[NUM_PARTICLES];

	glm::vec2 pos_offset(float(DOMAIN_WIDTH) * 0.5f, float(DOMAIN_HEIGHT) * 0.5f);

	for (int i = 0; i < NUM_PARTICLES; i++)
	{
		particle_positions[i] = glm::diskRand((float(DOMAIN_WIDTH) * 0.5f)) + pos_offset;
		particle_indices[i] = i;
	}

	glBindVertexArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, particles_indices_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, buff_size, (const GLvoid*)particle_indices);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_ARRAY_BUFFER, particles_position_vbo);

	buff_size = 0;

	for (GLint i = 0; i < NUM_POSITIONS_STORED; i++)
	{
		glBufferSubData(GL_ARRAY_BUFFER, buff_size, BUFFER_UPDATE_SIZE, (const GLvoid*)glm::value_ptr(particle_positions[0]));
		buff_size += BUFFER_UPDATE_SIZE;
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	delete[] particle_positions;
	delete[] particle_indices;
	particle_positions = 0;
	particle_indices = 0;
}

void DestroyParticleBuffers(void)
{
	glBindVertexArray(0);
	if (particles_position_vbo != 0)
	{
		glDeleteBuffers(1, &particles_position_vbo);
		particles_position_vbo = 0;
	}

	if (particles_indices_vbo != 0)
	{
		glDeleteBuffers(1, &particles_indices_vbo);
		particles_indices_vbo = 0;
	}

	if (particles_dest_vbo != 0)
	{
		glDeleteBuffers(1, &particles_dest_vbo);
		particles_dest_vbo = 0;
	}
}

void CreateParticleUpdateArrays(void)
{
	glBindVertexArray(0);

	if (!glIsBuffer(particles_position_vbo))
	{
		CreateParticleBuffers();
		FillParticleBuffers();
	}

	if (particles_update_vao != 0)
	{
		glDeleteVertexArrays(1, &particles_update_vao);
		particles_update_vao = 0;
	}

	glGenVertexArrays(1, &particles_update_vao);

	glBindVertexArray(particles_update_vao);

	glBindBuffer(GL_ARRAY_BUFFER, particles_position_vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);

	glBindVertexArray(0);
}

void DestroyParticleUpdateArrays(void)
{
	glBindVertexArray(0);

	if (particles_update_vao != 0)
	{
		glDeleteVertexArrays(1, &particles_update_vao);
		particles_update_vao = 0;
	}

}

void CreateVelocityTexture(void)
{
	if (!glIsBuffer(grid_vel_vbo))
	{
		printf("Error: Velocity Buffer Needs to be initalized before creating velocity texture\n");
		return;
	}

	if (glIsTexture(velocity_tex))
	{
		DestroyVelocityTexture();
	}

	glActiveTexture(GL_TEXTURE0);

	glGenTextures(1, &velocity_tex);

	glBindTexture(GL_TEXTURE_BUFFER, velocity_tex);

	glTexBuffer(GL_TEXTURE_BUFFER, GL_RG32F, grid_vel_vbo);

	glBindTexture(GL_TEXTURE_BUFFER, 0);
}

void DestroyVelocityTexture(void)
{
	if (glIsTexture(velocity_tex))
	{
		glDeleteTextures(1, &velocity_tex);
		velocity_tex = 0;
	}
}

void CreatePositionTexture(void)
{
	if (!glIsBuffer(particles_position_vbo))
	{
		printf("Error: Position Buffer Needs to be initalized before creating velocity texture\n");
		return;
	}

	if (glIsTexture(particles_position_tex))
	{
		DestroyPositionTexture();
	}

	glActiveTexture(GL_TEXTURE0);

	glGenTextures(1, &particles_position_tex);

	glBindTexture(GL_TEXTURE_BUFFER, particles_position_tex);

	glTexBuffer(GL_TEXTURE_BUFFER, GL_RG32F, particles_position_vbo);

	glBindTexture(GL_TEXTURE_BUFFER, 0);
}

void DestroyPositionTexture(void)
{
	if (glIsTexture(particles_position_tex))
	{
		glDeleteTextures(1, &particles_position_tex);
		particles_position_tex = 0;
	}
}

void Init(int argc, char** argv)
{
	GLenum glewInitResult;
	float lineWidthGran[2];
	float lineWidthRange[2];
	float pointFadeThresh[1];
	float pointSizeGran[1];
	int numSampleBuffers[1];
	GLint glutVersion, maxDrawBuffs, maxTexSize, maxColorAttachments, majVersion, minVersion, num_extensions, maxUniformBlockBindings, maxUniformBlockSize;

	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(CurrentWidth, CurrentHeight);

	MainWindow = glutCreateWindow(WINDOW_TITLE_PREFIX);

	glutVersion = glutGet(GLUT_VERSION);

	glewInitResult = glewInit();

	if (GLEW_OK != glewInitResult)
	{
		fprintf(stderr, "ERROR: %s\n", glewGetErrorString(glewInitResult));
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "INFO: GLUT Version: %d\n", glutVersion);
	fprintf(stdout, "INFO: OpenGL Version: %s\n", glGetString(GL_VERSION));
	fprintf(stdout, "INFO: OpenGL Renderer: %s\n", glGetString(GL_RENDERER));
	fprintf(stdout, "INFO: OpenGL Vendor: %s\n", glGetString(GL_VENDOR));
	fprintf(stdout, "INFO: OpenGL Shading Language Version: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

	glGetIntegerv(GL_MAJOR_VERSION, &majVersion);
	glGetIntegerv(GL_MINOR_VERSION, &minVersion);
	glGetIntegerv(GL_MAX_DRAW_BUFFERS, &maxDrawBuffs);
	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxTexSize);
	glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS, &maxColorAttachments);
	glGetIntegerv(GL_NUM_EXTENSIONS, &num_extensions);
	glGetIntegerv(GL_MAX_UNIFORM_BUFFER_BINDINGS, &maxUniformBlockBindings);
	glGetIntegerv(GL_MAX_UNIFORM_BLOCK_SIZE, &maxUniformBlockSize);

	glGetIntegerv(GL_SAMPLES, numSampleBuffers);
	glGetFloatv(GL_LINE_WIDTH_GRANULARITY, lineWidthGran);
	glGetFloatv(GL_LINE_WIDTH_RANGE, lineWidthRange);
	glGetFloatv(GL_POINT_FADE_THRESHOLD_SIZE, pointFadeThresh);
	glGetFloatv(GL_POINT_SIZE_GRANULARITY, pointSizeGran);

	glutInitContextVersion(majVersion, minVersion);
	glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);

	fprintf(stdout, "OPENGL VERSION: %d %d\n", majVersion, minVersion);
	fprintf(stdout, "NUM EXTENSIONS: %d\n", num_extensions);
	fprintf(stdout, "MAX DRAW BUFFERS SUPPORTED: %d\n", maxDrawBuffs);
	fprintf(stdout, "MAX TEXTURE SIZE SUPPORTED: %d\n", maxTexSize);
	fprintf(stdout, "MAX COLOR ATTACHMENTS SUPPORTED: %d\n", maxColorAttachments);
	fprintf(stdout, "MAX UNIFORM BUFFER BINDINGS: %d\n", maxUniformBlockBindings);
	fprintf(stdout, "MAX UNIFORM BLOCK SIZE: %d\n", maxUniformBlockSize);
	fprintf(stdout, "GL_SAMPLES: %d\n", numSampleBuffers[0]);
	fprintf(stdout, "GL_LINE_WIDTH_GRANULARITY: %f, %f\n", lineWidthGran[0], lineWidthGran[1]);
	fprintf(stdout, "GL_LINE_WIDTH_RANGE: %f, %f\n", lineWidthRange[0], lineWidthRange[1]);
	fprintf(stdout, "GL_POINT_FADE_THRESHOLD_SIZE: %f\n", pointFadeThresh[0]);
	fprintf(stdout, "GL_POINT_SIZE_GRANULARITY: %f\n", pointSizeGran[0]);

	fprintf(stdout, "INFO: OpenGL Version: %s\n", glGetString(GL_VERSION));

	glutDisplayFunc(Render);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);
	glutSpecialFunc(SpecialKey);
	glutTimerFunc(0, Timer, 0);
	glutIdleFunc(Idle);
	glutMouseFunc(Mouse);
	glutMotionFunc(MouseMove);

	glViewport(0, 0, CurrentWidth, CurrentHeight);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	printf("Attempting to enable GL_POINT_SPRITE\n");
	glEnable(GL_POINT_SPRITE);
	printf("Attempting to enable GL_BLEND\n");
	glEnable(GL_BLEND);
	printf("Attempting to enable GL_LINE_SMOOTH\n");
	glEnable(GL_LINE_SMOOTH);
	printf("Attempting to enable GL_LINE_SMOOTH_HINT\n");
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	printf("Attempting to enable GL_DEPTH_FUNC\n");
	glDepthFunc(GL_LEQUAL);
	printf("Attempting to enable GL_BLEND_FUNC\n");
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	printf("Attempting to enable GL_POINT_PARAMETER\n");
	glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
}

void UpdateParticles(void)
{
	GLuint startIndex = (GLuint)NUM_PARTICLES * CURRENT_INDEX_OFFSET;

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_BUFFER, velocity_tex);

	glUseProgram(particleUpdatesProgHandle);

	glUniformMatrix3fv(update_trans_loc, 1, GL_FALSE, glm::value_ptr(UpdateTransform));
	glUniformMatrix3fv(update_inv_trans_loc, 1, GL_FALSE, glm::value_ptr(UpdateInverseTransform));
	glUniform2fv(update_vel_scale_loc, 1, glm::value_ptr(velocity_scales));
	glUniform1f(update_dt_loc, dt);
	glUniform1i(update_width_loc, (VEL_GRID_RES_X + 1));
	glUniform1i(update_height_loc, (VEL_GRID_RES_Y + 1));
	glUniform1i(update_tex_loc, 0);

	glBindVertexArray(particles_update_vao);

	glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, particles_dest_vbo);

	glEnable(GL_RASTERIZER_DISCARD);

	glBeginTransformFeedback(GL_POINTS);
	glDrawArrays(GL_POINTS, startIndex, NUM_PARTICLES);
	glEndTransformFeedback();

	glDisable(GL_RASTERIZER_DISCARD);

	glBindVertexArray(0);

	glUseProgram(0);

	glBindTexture(GL_TEXTURE_BUFFER, 0);

	CURRENT_INDEX_OFFSET = (CURRENT_INDEX_OFFSET + 1) % NUM_POSITIONS_STORED;

	CopyParticleBuffers();
}

void CopyParticleBuffers(void)
{
	CURRENT_STORE_OFFSET = (CURRENT_INDEX_OFFSET * BUFFER_UPDATE_SIZE);

	glBindBuffer(GL_COPY_WRITE_BUFFER, particles_position_vbo);
	glBindBuffer(GL_COPY_READ_BUFFER, particles_dest_vbo);

	glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, CURRENT_STORE_OFFSET, BUFFER_UPDATE_SIZE);

	glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
	glBindBuffer(GL_COPY_READ_BUFFER, 0);
}

void InitUniforms(void)
{
	float left = 0.0f;
	float bottom = 0.0f;
	float right = float(CurrentWidth);
	float top = float(CurrentHeight);
	float scaleWidth, scaleHeight;
	vec3 col0, col1, col2;

	RenderMVP = glm::ortho(left, right, bottom, top);

	left = float(DOMAIN_MIN_X);
	right = float(DOMAIN_MAX_X);
	bottom = float(DOMAIN_MIN_Y);
	top = float(DOMAIN_MAX_Y);

	ParticlesMVP = glm::ortho(left, right, bottom, top);

	velocity_scales = vec2(0.05f, 0.05f);

	scaleWidth = right - left;
	scaleHeight = top - bottom;

	col0 = vec3(1.0f / scaleWidth, 0.0f, 0.0f);
	col1 = vec3(0.0f, 1.0f / scaleHeight, 0.0f);
	col2 = vec3(-left / scaleWidth, -bottom / scaleHeight, 1.0f);

	UpdateTransform = mat3(col0, col1, col2);

	col0 = vec3(scaleWidth, 0.0f, 0.0f);
	col1 = vec3(0.0f, scaleHeight, 0.0f);
	col2 = vec3(left, bottom, 1.0f);

	UpdateInverseTransform = mat3(col0, col1, col2);

	//	printf("RenderMVP:\n%s\n", glm::to_string(RenderMVP).c_str());
	//	printf("ParticlesMVP:\n%s\n", glm::to_string(ParticlesMVP).c_str());
	//	printf("Update Transform:\n%s\n", glm::to_string(UpdateTransform).c_str());
	//	printf("Update Inverse Transform:\n%s\n", glm::to_string(UpdateInverseTransform).c_str());
}

void Mouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)  // mouse pressed
	{
		if (force_mode == 1)
		{
			drag_path.clear();
			drag_path.push_back(glm::ivec2(x, (CurrentHeight - y)));
		}
		else
		{
			mouse_p1.x = double(x) / double(CurrentWidth);
			mouse_p1.y = double(CurrentHeight - y) / double(CurrentHeight);
		}

		left_mouse_down = true;
	}
	else // mouse released
	{
		if (force_mode == 1)
		{
			double invDX = 1.0 / double(CurrentWidth);
			double invDY = 1.0 / double(CurrentHeight);
			std::vector<dvec4> force_path;
			std::vector<ivec2>::const_iterator it;
			for (it = drag_path.begin(); it != drag_path.end(); ++it)
			{
				double x = double(it->x) * invDX;
				double y = double(it->y) * invDY;
				force_path.push_back(dvec4(x, y, 0.0, 0.0));
			}

			unsigned int sz = force_path.size() - 1;
			for (unsigned int idx = 0; idx < sz; idx++)
			{
				double dx = (force_path[idx + 1].x - force_path[idx].x) * force_mult;
				double dy = (force_path[idx + 1].y - force_path[idx].y) * force_mult;
				force_path[idx].z = dx;
				force_path[idx].w = dy;
			}

			Stir(force_path);
		}

		left_mouse_down = false;
	}
}

void MouseMove(int x, int y)
{
	if (left_mouse_down)
	{
		if (force_mode == 1)
		{
			drag_path.push_back(ivec2(x, (CurrentHeight - y)));
		}
		else
		{
			mouse_p2.x = double(x) / double(CurrentWidth);
			mouse_p2.y = double(CurrentHeight - y) / double(CurrentHeight);

			double fx = (mouse_p2.x - mouse_p1.x) * force_mult;
			double fy = (mouse_p2.y - mouse_p1.y) * force_mult;

			std::vector<dvec4> force_path(2, dvec4(0.0));
			force_path[0] = glm::dvec4(mouse_p1.x, mouse_p1.y, fx, fy);
			Stir(force_path);
			mouse_p1 = mouse_p2;
		}
	}
}

void Keyboard(unsigned char key, int x, int y)
{
	int n;
	ivec2 k;
	switch (key)
	{
	case 27:
		glutLeaveMainLoop();
		break;
	case '=':
		viscosity *= 1.01f;
		printf("viscosity: %f\n", viscosity);
		break;
	case '-':
		viscosity *= 0.99f;
		printf("viscosity: %f\n", viscosity);
		break;
	case ']':
		dt *= 1.01f;
		printf("dt: %f\n", dt);
		break;
	case '[':
		dt *= 0.99f;
		printf("dt: %f\n", dt);
		break;
	case 'r':
		printf("Resetting\n");
		Reset();
		break;
	case 'b':
		n = (CURRENT_RENDER_MODE + 1) % RENDER_NUM_MODES;
		SetRenderMode(n);
		switch (CURRENT_RENDER_MODE)
		{
		case RENDER_VELOCITY_GRID:
			printf("Render Mode: Velocity Grid\n");
			break;
		case RENDER_PARTICLES:
			printf("Render Mode: Particles\n");
			break;
		case RENDER_GRID_AND_PARTICLES:
			printf("Render Mode: Grid And Particles\n");
			break;
		default:
			printf("Render Mode: NONE\n");
			break;
		}
		break;
	case 'v':
		UPDATE_VELOCITY_GRID_ON = !UPDATE_VELOCITY_GRID_ON;
		printf("UPDATE_VELOCITY_GRID_ON: %d\n", UPDATE_VELOCITY_GRID_ON);
		break;
	case 'e':
		UPDATE_PARTICLES_ON = !UPDATE_PARTICLES_ON;
		printf("UPDATE_PARTICLES_ON: %d\n", UPDATE_PARTICLES_ON);
		break;
	case 'c':
		printf("Updating Particles: CURRENT_INDEX_OFFSET: %d CURRENT_SIZE_OFFSET: %d\n", CURRENT_INDEX_OFFSET, CURRENT_STORE_OFFSET);
		UpdateParticles();
		printf("CURRENT_INDEX_OFFSET: %d CURRENT_SIZE_OFFSET: %d\n", CURRENT_INDEX_OFFSET, CURRENT_STORE_OFFSET);
		break;
	case '.':
		force_mult *= 1.01f;
		printf("force_mult: %f\n", force_mult);
		break;
	case ',':
		force_mult *= 0.99f;
		printf("force_mult: %f\n", force_mult);
		break;
	case 'l':
		n = (N_SQ + 1) * (N_SQ + 1);
		SetBasisDimension(n);
		break;
	case 'k':
		n = (N_SQ - 1) * (N_SQ - 1);
		SetBasisDimension(n);
		break;
	case 'f':
		if (force_mode == 2)
		{
			force_mode = 1;
		}
		else
		{
			force_mode = 2;
		}
		printf("Force Mode: %d\n", force_mode);
		break;
	case 's':
		n = int(glm::floor(sqrt(double(VEL_GRID_RES_X)))) + 1;
		SetVelocityResolution((n * n));
		printf("Velocity Resolution %d\n", n);
		break;
	case 'a':
		n = int(glm::floor(sqrt(double(VEL_GRID_RES_Y)))) - 1;
		SetVelocityResolution((n * n));
		printf("Velocity Resolution %d\n", n);
		break;
	case 'z':
		if (viscosity == 0.0)
		{
			viscosity = 0.001;
		}
		else
		{
			viscosity = 0.0;
		}
		printf("viscosity: %f\n", viscosity);
		break;
	case 'w':
		ZERO_MODES_ON = !ZERO_MODES_ON;
		printf("zero modes: %d\n", ZERO_MODES_ON);
		break;
	case 'p':
		k = BasisLookupK(CURRENT_MODE);
		printf("Current Mode %d = [%d,%d]\n", CURRENT_MODE, k.x, k.y);
		break;
	case 'm':
		particle_tail_line_width += .25f;
		if (particle_tail_line_width > 63.0f)
		{
			particle_tail_line_width = 63.0f;
		}
		printf("particle_tail_line_width: %f\n", particle_tail_line_width);
		break;
	case 'n':
		particle_tail_line_width -= 0.25f;
		if (particle_tail_line_width < 0.25f)
		{
			particle_tail_line_width = 0.25f;
		}
		printf("particle_tail_line_width: %f\n", particle_tail_line_width);
		break;
	}
}

void SpecialKey(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_LEFT:
		if (CURRENT_MODE == 0)
		{
			CURRENT_MODE = N - 1;
		}
		else
		{
			CURRENT_MODE = (CURRENT_MODE - 1) % (GLuint)N;
		}
		printf("CurrentMode: %d\n", CURRENT_MODE);
		break;
	case GLUT_KEY_RIGHT:
		CURRENT_MODE = (CURRENT_MODE + 1) % (GLuint)N;
		printf("CurrentMode: %d\n", CURRENT_MODE);
		break;
	case GLUT_KEY_UP:
		printf("Inserting Mode: %d\n", CURRENT_MODE);
		zero_modes[CURRENT_MODE] = 1;
		break;
	case GLUT_KEY_DOWN:
		printf("Zeroing Mode: %d\n", CURRENT_MODE);
		zero_modes[CURRENT_MODE] = 0;
		break;
	case GLUT_KEY_HOME:
		printf("Inserting Multiples Of: %d\n", CURRENT_MODE);
		InsertMultiplesOf(CURRENT_MODE);
		break;
	case GLUT_KEY_END:
		printf("Zeroing Multiples Of: %d\n", CURRENT_MODE);
		ZeroMultiplesOf(CURRENT_MODE);
		break;
	}
}

void Timer(int Value)
{
	if (0 != Value)
	{
		//	char* TempString = (char*)malloc(512 + strlen(WINDOW_TITLE_PREFIX));

		//	sprintf_s(TempString, 1, "%s: %d Frames Per Second @ %d x %d", WINDOW_TITLE_PREFIX, FrameCount * 4, CurrentWidth, CurrentHeight);
		////	sprintf_s(TempString, "%s: %d Frames Per Second @ %d x %d", WINDOW_TITLE_PREFIX, FrameCount * 4, CurrentWidth, CurrentHeight);

		//	glutSetWindowTitle(TempString);
		//	free(TempString);

		char* TempString = (char*)malloc(512 + strlen(WINDOW_TITLE_PREFIX));
		// Correct the buffer size to match the allocated size
		sprintf_s(TempString, 512 + strlen(WINDOW_TITLE_PREFIX), "%s: %d Frames Per Second @ %d x %d", WINDOW_TITLE_PREFIX, FrameCount * 4, CurrentWidth, CurrentHeight);

		glutSetWindowTitle(TempString);
		free(TempString);
	}
	FrameCount = 0;
	glutTimerFunc(250, Timer, 1);
}

void Idle(void)
{
	if (vel_buff_size_changed)
	{
		grid_pos_needs_resize = true;
		grid_vel_needs_resize = true;
	}

	if (ZERO_MODES_ON)
	{
		ZeroModes();
	}

	if (UPDATE_VELOCITY_GRID_ON)
	{
		UpdateFields();
	}

	if (UPDATE_PARTICLES_ON)
	{
		UpdateParticles();
	}

	if (vel_buff_size_changed)
	{
		if ((!grid_vel_needs_resize) && (!grid_pos_needs_resize))
		{
			vel_buff_size_changed = false;
		}
	}

	glutPostRedisplay();
}

void Reshape(int w, int h)
{
	CurrentWidth = w;
	CurrentHeight = h;

	float left = 0.0f;
	float down = 0.0f;
	float right = float(CurrentWidth);
	float up = float(CurrentHeight);

	RenderMVP = glm::ortho(left, right, down, up);

	glutPostRedisplay();
}

void InitFields(void)
{
	wCoefficients.set_size(N);
	wCoefficients.zeros();
	dwCoefficients.set_size(N);
	dwCoefficients.zeros();
	dwForceCoefficients.set_size(N);
	dwForceCoefficients.zeros();

	for (int i = 0; i < 4; i++)
	{
		rk4Qn[i].set_size(N);
		rk4Qn[i].zeros();
		rk4Dwt[i].set_size(N);
		rk4Dwt[i].zeros();
	}

	wCoefficients(0) = 1.0;

	CreateCkMatrices();
}

void FillLookupTable(void)
{
	CreateLookupTables();

	int index = 0;
	for (int k1 = 0; k1 <= N_SQ; k1++)
	{
		for (int k2 = 0; k2 <= N_SQ; k2++)
		{
			if (k1 > N_SQ || k1 < 1 || k2 > N_SQ || k2 < 1)
			{
				continue;
			}

			basis_lookup_table[index].x = k1;
			basis_lookup_table[index].y = k2;

			basis_rlookup_table[k1][k2] = index;
			index++;
		}
	}
}

void PrecomputeBasisFields(void)
{
	CreateBasisField();

	for (int k = 0; k < N; k++)
	{
		glm::ivec2 idxK;
		idxK.x = BasisLookup(k, 0);
		idxK.y = BasisLookup(k, 1);
		BasisFieldRect2D(idxK, initial_amp_scale, velocity_basis[k]);
	}
}

void PrecomputeDynamics(void)
{
	CreateEigenValues();

	for (int k = 0; k < N; k++)
	{
		glm::ivec2 idxK;
		idxK.x = BasisLookup(k, 0);
		idxK.y = BasisLookup(k, 1);

		double eVal = double(idxK.x * idxK.x + idxK.y * idxK.y);

		eigs[k] = eVal;
		eigs_inv[k] = 1.0 / eVal;
		eigs_inv_root[k] = 1.0 / sqrt(eVal);
	}

	for (int d1 = 0; d1 < N; d1++)
	{
		glm::ivec2 idxA;
		idxA.x = BasisLookup(d1, 0);
		idxA.y = BasisLookup(d1, 1);

		double lambda_a = -double(idxA.x * idxA.x + idxA.y * idxA.y);

		for (int d2 = 0; d2 < N; d2++)
		{
			glm::ivec2 idxB;
			idxB.x = BasisLookup(d2, 0);
			idxB.y = BasisLookup(d2, 1);

			double lambda_b = -double(idxB.x * idxB.x + idxB.y * idxB.y);
			double inv_lambda_b = 1.0 / lambda_b;

			int k1 = BasisRLookup(idxA);
			int k2 = BasisRLookup(idxB);

			glm::ivec2 antipairs[4];

			antipairs[0].x = idxA.x - idxB.x;
			antipairs[1].x = idxA.x - idxB.x;
			antipairs[2].x = idxA.x + idxB.x;
			antipairs[3].x = idxA.x + idxB.x;

			antipairs[0].y = idxA.y - idxB.y;
			antipairs[1].y = idxA.y + idxB.y;
			antipairs[2].y = idxA.y - idxB.y;
			antipairs[3].y = idxA.y + idxB.y;

			for (int c = 0; c < 4; c++)
			{
				int index = BasisRLookup(antipairs[c]);

				if (index != -1)
				{
					glm::ivec2 idxC(c, 0);
					double coef = -CoefDensity(idxA, idxB, idxC) * inv_lambda_b;
					double coef2 = -(lambda_b / lambda_a) * coef;

					if (coef != 0.0)
					{
#ifdef WIN32
						CkMatrices[index](k1, k2) = coef;
#else
						CkMatrices[index].insert(k1, k2) = coef;
#endif
					}
					if (coef2 != 0.0)
					{
#ifdef WIN32
						CkMatrices[index](k2, k1) = coef2;
#else
						CkMatrices[index].insert(k2, k1) = coef2;
#endif
					}

				}
			}
		}
	}
}

void BasisFieldRect2D(const glm::ivec2& K, double amp, glm::dvec2** vBasisElem)
{
	double scaleX = 1.0;
	double scaleY = 1.0;

	if (K.x != 0)
	{
		scaleX = -1.0 / double(K.x * K.x + K.y * K.y);
	}
	if (K.y != 0)
	{
		scaleY = -1.0 / double(K.x * K.x + K.y * K.y);
	}

	double dx = DOMAIN_WIDTH / double(VEL_GRID_RES_X);
	double dy = DOMAIN_HEIGHT / double(VEL_GRID_RES_Y);
	double k1 = double(K.x);
	double k2 = double(K.y);

	for (int row = 0; row <= VEL_GRID_RES_Y; row++)
	{
		double y0 = double(row) * dy;
		double y1 = y0 + dy * 0.5;
		double sinY0K2 = sin(y0 * k2);
		double cosY1K2 = cos(y1 * k2);

		for (int col = 0; col <= VEL_GRID_RES_X; col++)
		{
			double x0 = double(col) * dx;
			double x1 = x0 + dx * 0.5;
			double sinX0K1 = sin(x0 * k1);
			double cosX1K1 = cos(x1 * k1);

			vBasisElem[row][col].x = -amp * scaleX * k2 * sinX0K1 * cosY1K2;
			vBasisElem[row][col].y = amp * scaleY * k1 * cosX1K1 * sinY0K2;
		}
	}
}

void ExpandBasis(void)
{
	GLfloat* ptr = 0;
	glBindBuffer(GL_ARRAY_BUFFER, grid_vel_vbo);

	ptr = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

	for (int row = 0; row <= VEL_GRID_RES_Y; row++)
	{
		for (int col = 0; col <= VEL_GRID_RES_X; col++)
		{
			dvec2 vSum = dvec2(0.0, 0.0);
			for (int k = 0; k < N; k++)
			{
				vSum += (velocity_basis[k][row][col] * wCoefficients(k));
			}

			*ptr = float(vSum.x);
			ptr++;
			*ptr = float(vSum.y);
			ptr++;
		}
	}

	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ProjectForces(const std::vector<dvec4>& force_path, arma::vec& delW)
{
	delW.set_size(N);
	delW.zeros();
	for (int k = 0; k < N; k++)
	{
		int k1 = BasisLookup(k, 0);
		int k2 = BasisLookup(k, 1);

		double xfact = 1.0;
		double yfact = 1.0;
		if (k1 != 0)
		{
			xfact = -1.0 / (double(k1 * k1 + k2 * k2));
		}
		if (k2 != 0)
		{
			yfact = -1.0 / (double(k1 * k1 + k2 * k2));
		}

		std::vector<dvec4>::const_iterator it;
		for (it = force_path.begin(); it != force_path.end(); ++it)
		{
			double x = it->x;
			double y = it->y;
			double fx = it->z;
			double fy = it->w;

			if (x >= 1.00001 || x <= -0.00001 || y >= 1.00001 || y <= -0.00001)
			{
				continue;
			}

			x *= PI;
			y *= PI;

			double vx = -double(k2) * xfact * sin(double(k1) * x) * cos(double(k2) * y) * dt;
			double vy = double(k1) * yfact * cos(double(k1) * x) * sin(double(k2) * y) * dt;

			delW(k) += (vx * fx + vy * fy);
		}
	}
}

void Stir(const std::vector<dvec4>& force_path)
{
	arma::vec delW;
	ProjectForces(force_path, delW);
	dwForceCoefficients += delW;
}

double CurrentEnergy(void)
{
	double energy = 0.0;
	for (int k = 0; k < N; k++)
	{
		double c2 = wCoefficients(k) * wCoefficients(k);
		energy += (eigs_inv[k] * c2);
	}
	return energy;
}

void SetEnergy(double desired_e)
{
	double cur_e = CurrentEnergy();
	double fact = sqrt(desired_e) / sqrt(cur_e);

	wCoefficients = wCoefficients * fact;
}

int BasisLookup(int index, int component)
{
	return basis_lookup_table[index][component];
}

ivec2 BasisLookupK(int index)
{
	return basis_lookup_table[index];
}

int BasisRLookup(const glm::ivec2& K)
{
	if (K.x > N_SQ || K.x < 1 || K.y > N_SQ || K.y < 1)
	{
		return -1;
	}
	return basis_rlookup_table[K.x][K.y];
}

double CoefDensity(const glm::ivec2& a, const glm::ivec2& b, const glm::ivec2& c)
{
	switch (c.y)
	{
	case 0:
		switch (c.x)
		{
		case 0:
			return -double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return -double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		}
		break;
	case 1:
		switch (c.x)
		{
		case 0:
			return -double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return -double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		}
		break;
	case 2:
		switch (c.x)
		{
		case 0:
			return -double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return -double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		}
		break;
	case 3:
		switch (c.x)
		{
		case 0:
			return -double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return -double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return double(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return double(a.x * b.y - b.x * a.y) * 0.25;
			break;
		}
		break;
	}

	return 0;
}

void CreateLookupTables(void)
{
	if (basis_lookup_table != NULL)
	{
		DestroyLookupTables();
	}

	N_SQ = int(glm::floor(sqrt(double(N))));

	basis_lookup_table = new glm::ivec2[N];
	basis_rlookup_table = new int* [N_SQ + 1];
	for (int i = 0; i < N; i++)
	{
		basis_lookup_table[i].x = -1;
		basis_lookup_table[i].y = -1;
	}

	for (int i = 0; i <= N_SQ; i++)
	{
		basis_rlookup_table[i] = new int[N_SQ + 1];
	}

	for (int k1 = 0; k1 <= N_SQ; k1++)
	{
		for (int k2 = 0; k2 <= N_SQ; k2++)
		{
			basis_rlookup_table[k1][k2] = -1;
		}
	}
}

void CreateCkMatrices(void)
{
	if (CkMatrices != NULL)
	{
		DestroyCkMatrices();
	}

	CkMatrices = new arma::sp_mat[N];

	for (int k = 0; k < N; k++)
	{
		CkMatrices[k].set_size(N, N);
	}
}

void CreateBasisField(void)
{
	if (velocity_basis != NULL)
	{
		DestroyBasisField();
	}

	int* tmp_modes = new int[N];

	for (int i = 0; i < N; i++)
	{
		if (i < N_OLD && zero_modes != NULL)
		{
			tmp_modes[i] = zero_modes[i];
		}
		else
		{
			tmp_modes[i] = 1;
		}
	}

	if (zero_modes != NULL)
	{
		delete[] zero_modes;
		zero_modes = NULL;
	}

	zero_modes = tmp_modes;

	velocity_basis = new glm::dvec2 * *[N];

	for (int k = 0; k < N; k++)
	{
		velocity_basis[k] = new glm::dvec2 * [VEL_GRID_RES_Y + 1];
		for (int row = 0; row <= VEL_GRID_RES_Y; row++)
		{
			velocity_basis[k][row] = new glm::dvec2[VEL_GRID_RES_X + 1];
			for (int col = 0; col <= VEL_GRID_RES_X; col++)
			{
				velocity_basis[k][row][col] = glm::dvec2(0.0, 0.0);
			}
		}
	}
}

void CreateEigenValues(void)
{
	if (eigs != NULL)
	{
		DestroyEigenValues();
	}

	eigs = new double[N];
	eigs_inv = new double[N];
	eigs_inv_root = new double[N];

	for (int k = 0; k < N; k++)
	{
		eigs[k] = 0.0;
		eigs_inv[k] = 0.0;
		eigs_inv_root[k] = 0.0;
	}
}

void DestroyLookupTables(void)
{
	if (basis_lookup_table != NULL)
	{
		delete[] basis_lookup_table;
		basis_lookup_table = NULL;
	}

	if (basis_rlookup_table != NULL)
	{
		for (int i = 0; i <= N_SQ; i++)
		{
			if (basis_rlookup_table[i] != NULL)
			{
				delete[] basis_rlookup_table[i];
				basis_rlookup_table[i] = NULL;
			}
		}

		delete[] basis_rlookup_table;
		basis_rlookup_table = NULL;
	}
}

void DestroyCkMatrices(void)
{
	if (CkMatrices != NULL)
	{
		delete[] CkMatrices;
		CkMatrices = NULL;
	}
}

void DestroyBasisField(void)
{
	if (velocity_basis != NULL)
	{
		for (int k = 0; k < N; k++)
		{
			if (velocity_basis[k] != NULL)
			{
				for (int row = 0; row <= VEL_GRID_RES_Y; row++)
				{
					if (velocity_basis[k][row] != NULL)
					{
						delete[] velocity_basis[k][row];
						velocity_basis[k][row] = NULL;
					}
				}

				delete[] velocity_basis[k];
				velocity_basis[k] = NULL;
			}
		}

		delete[] velocity_basis;
		velocity_basis = NULL;
	}
}

void DestroyEigenValues(void)
{
	if (eigs != NULL)
	{
		delete[] eigs;
		eigs = NULL;
	}
	if (eigs_inv != NULL)
	{
		delete[] eigs_inv;
		eigs_inv = NULL;
	}
	if (eigs_inv_root != NULL)
	{
		delete[] eigs_inv_root;
		eigs_inv_root = NULL;
	}
}

void Shutdown(void)
{
	DestroyGrid();
	DestroyParticleRenderArrays();
	DestroyParticleBuffers();
	DestroyVelocityTexture();
	DestroyParticleUpdateArrays();
	DestroyLookupTables();
	DestroyBasisField();
	DestroyEigenValues();
	DestroyCkMatrices();
}

void SetVelocityResolution(int res)
{
	if (res <= 0)
	{
		printf("Error: Cannot set velocity resolution %d to non positive integer\n", res);
		return;
	}

	printf("Attempting to set velocity resolution to: %d\n", res);

	DestroyLookupTables();
	DestroyCkMatrices();
	DestroyBasisField();
	DestroyEigenValues();

	VEL_GRID_RES_X = res;
	VEL_GRID_RES_Y = res;

	vel_buff_size_changed = true;
	grid_pos_needs_resize = true;
	grid_vel_needs_resize = true;

	InitFields();
	FillLookupTable();
	PrecomputeBasisFields();
	PrecomputeDynamics();
	ExpandBasis();

	FillGridPositionBuffer();
	FillGridVelocityBuffer();
	DestroyVelocityTexture();
	CreateVelocityTexture();
	FillParticleBuffers();
}

void SetBasisDimension(int dim)
{
	if (dim <= 0)
	{
		printf("Error: Cannot set basis dimension to non positive integer\n");
		return;
	}

	printf("Attempting to set basis dimension to: %d\n", dim);

	DestroyLookupTables();
	DestroyCkMatrices();
	DestroyBasisField();
	DestroyEigenValues();

	N_OLD = N;

	N = dim;

	InitFields();
	FillLookupTable();
	PrecomputeBasisFields();
	PrecomputeDynamics();
	ExpandBasis();

	FillGridPositionBuffer();
	FillGridVelocityBuffer();
	DestroyVelocityTexture();
	CreateVelocityTexture();
	FillParticleBuffers();
}

void InitShaders(void)
{
	GLuint vs, gs, fs;
	GLint status;

	velocitiesProgHandle = glCreateProgram();
	vs = LoadShaderFromFile(GL_VERTEX_SHADER, "VelocityRender2D.vs.glsl");
	gs = LoadShaderFromFile(GL_GEOMETRY_SHADER, "VelocityRender2D.gs.glsl");
	fs = LoadShaderFromFile(GL_FRAGMENT_SHADER, "VelocityRender2D.fs.glsl");
	glAttachShader(velocitiesProgHandle, vs);
	glAttachShader(velocitiesProgHandle, gs);
	glAttachShader(velocitiesProgHandle, fs);
	glLinkProgram(velocitiesProgHandle);
	glGetProgramiv(velocitiesProgHandle, GL_LINK_STATUS, &status);
	if (status == GL_FALSE)
	{
		GLint infoLogLength;
		glGetProgramiv(velocitiesProgHandle, GL_INFO_LOG_LENGTH, &infoLogLength);

		GLchar* strInfoLog = new GLchar[infoLogLength + 1];
		glGetProgramInfoLog(velocitiesProgHandle, infoLogLength, NULL, strInfoLog);
		fprintf(stderr, "Linker failure: %s\n", strInfoLog);
		delete[] strInfoLog;
	}

	glDeleteShader(vs);
	glDeleteShader(gs);
	glDeleteShader(fs);

	particleTailsProgHandle = glCreateProgram();
	vs = LoadShaderFromFile(GL_VERTEX_SHADER, "ParticleTailsRender2D.vs.glsl");
	gs = LoadShaderFromFile(GL_GEOMETRY_SHADER, "ParticleTailsRender2D.gs.glsl");
	fs = LoadShaderFromFile(GL_FRAGMENT_SHADER, "ParticleTailsRender2D.fs.glsl");

	glAttachShader(particleTailsProgHandle, vs);
	glAttachShader(particleTailsProgHandle, gs);
	glAttachShader(particleTailsProgHandle, fs);
	glLinkProgram(particleTailsProgHandle);
	glGetProgramiv(particleTailsProgHandle, GL_LINK_STATUS, &status);
	if (status == GL_FALSE)
	{
		GLint infoLogLength;
		glGetProgramiv(particleTailsProgHandle, GL_INFO_LOG_LENGTH, &infoLogLength);

		GLchar* strInfoLog = new GLchar[infoLogLength + 1];
		glGetProgramInfoLog(particleTailsProgHandle, infoLogLength, NULL, strInfoLog);
		fprintf(stderr, "Linker failure: %s\n", strInfoLog);
		delete[] strInfoLog;
	}

	glDeleteShader(vs);
	glDeleteShader(gs);
	glDeleteShader(fs);

	const char* feedback_varyings[] = { { "NewPos2D" } };

	particleUpdatesProgHandle = glCreateProgram();
	vs = LoadShaderFromFile(GL_VERTEX_SHADER, "ParticleUpdate2D.vs.glsl");

	glAttachShader(particleUpdatesProgHandle, vs);

	glTransformFeedbackVaryings(particleUpdatesProgHandle, 1, feedback_varyings, GL_SEPARATE_ATTRIBS);

	glLinkProgram(particleUpdatesProgHandle);
	glGetProgramiv(particleUpdatesProgHandle, GL_LINK_STATUS, &status);
	if (status == GL_FALSE)
	{
		GLint infoLogLength;
		glGetProgramiv(particleUpdatesProgHandle, GL_INFO_LOG_LENGTH, &infoLogLength);

		GLchar* strInfoLog = new GLchar[infoLogLength + 1];
		glGetProgramInfoLog(particleUpdatesProgHandle, infoLogLength, NULL, strInfoLog);
		fprintf(stderr, "Linker failure: %s\n", strInfoLog);
		delete[] strInfoLog;
	}

	glDeleteShader(vs);

	vel_mvp_loc = glGetUniformLocation(velocitiesProgHandle, "MVP");
	vel_col_loc = glGetUniformLocation(velocitiesProgHandle, "LineColor");
	vel_len_scale_loc = glGetUniformLocation(velocitiesProgHandle, "VelocityLengthScale");

	particles_mvp_loc = glGetUniformLocation(particleTailsProgHandle, "MVP");
	particles_head_col_loc = glGetUniformLocation(particleTailsProgHandle, "HeadLineColor");
	particles_tail_col_loc = glGetUniformLocation(particleTailsProgHandle, "TailLineColor");
	particles_num_particles_loc = glGetUniformLocation(particleTailsProgHandle, "NumParticles");
	particles_num_stored_loc = glGetUniformLocation(particleTailsProgHandle, "NumPositionsStored");
	particles_current_store_loc = glGetUniformLocation(particleTailsProgHandle, "CurrentStoreOffset");
	particles_tex_buff_loc = glGetUniformLocation(particleTailsProgHandle, "StoredPositions2DBuffer");

	update_tex_loc = glGetUniformLocation(particleUpdatesProgHandle, "VelocityField");
	update_trans_loc = glGetUniformLocation(particleUpdatesProgHandle, "UpdateTransform");
	update_inv_trans_loc = glGetUniformLocation(particleUpdatesProgHandle, "UpdateInverseTransform");
	update_vel_scale_loc = glGetUniformLocation(particleUpdatesProgHandle, "VelocityScale");
	update_dt_loc = glGetUniformLocation(particleUpdatesProgHandle, "dt");
	update_width_loc = glGetUniformLocation(particleUpdatesProgHandle, "GridWidth");
	update_height_loc = glGetUniformLocation(particleUpdatesProgHandle, "GridHeight");

	printf("velocitiesProgHandle: %d mvp: %d col: %d len: %d\n", velocitiesProgHandle, vel_mvp_loc, vel_col_loc, vel_len_scale_loc);
	printf("particlesProgHandle: %d mvp: %d head col: %d tail col: %d num particles: %d num stored: %d current store: %d tex buff: %d\n", particleTailsProgHandle, particles_mvp_loc, particles_head_col_loc, particles_tail_col_loc, particles_num_particles_loc, particles_num_stored_loc, particles_current_store_loc, particles_tex_buff_loc);
	printf("updateProgHandle: %d tex: %d trans: %d inv_trans: %d vel: %d dt: %d width: %d height: %d\n", particleUpdatesProgHandle, update_tex_loc, update_trans_loc, update_inv_trans_loc, update_vel_scale_loc, update_dt_loc, update_width_loc, update_height_loc);
}

GLuint LoadShaderFromFile(GLenum type, const char* filename)
{
	assert((type == GL_VERTEX_SHADER) || (type == GL_FRAGMENT_SHADER) || (type == GL_GEOMETRY_SHADER) || (type == GL_TESS_CONTROL_SHADER) || (type == GL_TESS_EVALUATION_SHADER));
	GLuint shader_handle = 0;
	ifstream fp(filename, ios_base::in);
	if (fp)
	{
		string buffer(std::istreambuf_iterator<char>(fp), (std::istreambuf_iterator<char>()));
		shader_handle = glCreateShader(type);
		GLint shader_len = buffer.length();
		GLchar* tempShad = new GLchar[shader_len + 1];
		memcpy(tempShad, buffer.c_str(), shader_len);
		tempShad[shader_len] = '\0';
		glShaderSource(shader_handle, 1, (const GLchar**)&tempShad, NULL);
		delete[] tempShad;
		tempShad = 0;
		glCompileShader(shader_handle);
		//		printf("Loading Shader: %s\n", buffer.c_str());
		GLint compileStatus;
		glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compileStatus);
		if (compileStatus != GL_TRUE)
		{
			GLint infoLogLength;
			glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &infoLogLength);

			GLchar* strInfoLog = new GLchar[infoLogLength + 1];
			glGetShaderInfoLog(shader_handle, infoLogLength, NULL, strInfoLog);
			fprintf(stderr, "Compile Shader failure: %s\n", strInfoLog);
			delete[] strInfoLog;
		}
	}
	else
	{
		printf("Error loading shader: %s\n", filename);
	}
	return shader_handle;
}