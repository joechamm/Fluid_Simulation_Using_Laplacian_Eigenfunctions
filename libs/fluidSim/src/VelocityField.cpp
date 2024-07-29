//#include "VelocityField.h"
//#include "FluidSimulation.h"

//#include "../../include/VelocityField.h"
//#include "../../include/FluidSimulation.h"
#include <fluidSim/include/VelocityField.h>
#include <fluidSim/include/FluidSimulation.h>

#include <stdexcept>
#include <iostream>

VelocityField::VelocityField(uint32_t numGridCols, uint32_t numGridRows, uint32_t basisDimension, double amplitude)
	: m_basis_dimension(basisDimension),
	m_num_basis_functions(basisDimension * basisDimension),
	m_velocity_grid_resolution_x(numGridCols),
	m_velocity_grid_resolution_y(numGridRows),
	m_amplitude(amplitude),
	m_velocity_basis(nullptr)
{
	// ensure the basis dimension is at least 2
	try {
		if (m_basis_dimension < 2)
			throw std::invalid_argument("Basis dimension must be at least 2");

		// allocate memory for the velocity basis elements
		m_velocity_basis = new glm::dvec2**[m_num_basis_functions];
		for (uint32_t k = 0; k < m_num_basis_functions; k++)
		{
			m_velocity_basis[k] = new glm::dvec2 * [m_velocity_grid_resolution_y + 1];
			for (uint32_t row = 0; row <= m_velocity_grid_resolution_y; row++)
			{
				m_velocity_basis[k][row] = new glm::dvec2[m_velocity_grid_resolution_x + 1];
				for (uint32_t col = 0; col <= m_velocity_grid_resolution_x; col++)
				{
					m_velocity_basis[k][row][col] = glm::dvec2(0.0, 0.0);
				}
			}
		}
	}
	catch (const std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	}
}

VelocityField::~VelocityField()
{
	// deallocate memory for the velocity basis elements
	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		if (m_velocity_basis[k] != nullptr)
		{
			for (uint32_t row = 0; row <= m_velocity_grid_resolution_y; row++)
			{
				if (m_velocity_basis[k][row] != nullptr)
				{
					delete[] m_velocity_basis[k][row];
					m_velocity_basis[k][row] = nullptr;
				}
			}

			delete[] m_velocity_basis[k];
			m_velocity_basis[k] = nullptr;
		}		
	}
	delete[] m_velocity_basis;
	m_velocity_basis = nullptr;
}

void VelocityField::PrecomputeBasisField(const FluidSimulation& simulation)
{
	// make sure the simulation basis dimension matches the velocity field basis dimension
	try {
		if (simulation.GetBasisDimension() != m_basis_dimension)
			throw std::invalid_argument("Simulation basis dimension does not match velocity field basis dimension");

		for (uint32_t k = 0; k < m_num_basis_functions; k++)
		{
			glm::ivec2 idxK;
			idxK.x = simulation.BasisLookup(k, 0);
			idxK.y = simulation.BasisLookup(k, 1);
			simulation.BasisFieldRect2D(idxK, m_velocity_grid_resolution_x, m_velocity_grid_resolution_y, m_amplitude, m_velocity_basis[k]);
		}

	} catch(const std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	}
}

glm::dvec2 VelocityField::ComputeVelocity(uint32_t row, uint32_t col, const FluidSimulation& simulation) const
{
	glm::dvec2 vsum(0.0, 0.0);

	// make sure the simulation basis dimension matches the velocity field basis dimension
	try {
		for(uint32_t k = 0; k < m_num_basis_functions; k++)
		{
			vsum += m_velocity_basis[k][row][col] * simulation.GetWCoefficient(k);
		}
	}
	catch (const std::invalid_argument& e)
	{
		std::cerr << e.what() << std::endl;
	}

	return vsum;
}
