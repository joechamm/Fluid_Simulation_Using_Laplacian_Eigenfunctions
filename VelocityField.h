#pragma once

#include <glm/glm.hpp>

class FluidSimulation;

class VelocityField
{
	double m_amplitude;

	uint32_t m_basis_dimension;
	uint32_t m_num_basis_functions;

	uint32_t m_velocity_grid_resolution_x;
	uint32_t m_velocity_grid_resolution_y;
	glm::dvec2*** m_velocity_basis;

public:
	VelocityField(uint32_t numGridCols, uint32_t numGridRows, uint32_t basisDimension, double amplitude = 1.0);
	~VelocityField();

	/// Accessors
	uint32_t GetBasisDimension() const { return m_basis_dimension; }
	uint32_t GetNumBasisFunctions() const { return m_num_basis_functions; }
	uint32_t GetVelocityGridResolutionX() const { return m_velocity_grid_resolution_x; }
	uint32_t GetVelocityGridResolutionY() const { return m_velocity_grid_resolution_y; }

	double GetAmplitude() const { return m_amplitude; }
	void SetAmplitude(double amplitude) { m_amplitude = amplitude; }

	void PrecomputeBasisField(const FluidSimulation& simulation);
};

