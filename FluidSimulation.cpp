#include "FluidSimulation.h"

// TODO: set m_velocity_grid_resolution_y and m_velocity_grid_resolution_x 
FluidSimulation::FluidSimulation(uint32_t dimension, double viscosity, double dt, const glm::dvec4& domain) :
	m_basis_lookup_table(nullptr),
	m_basis_rlookup_table(nullptr),
	m_eigs(nullptr),
	m_eigs_inv(nullptr),
	m_eigs_inv_root(nullptr),
	m_CkMatrices(nullptr),
	m_zero_modes(nullptr),
	m_basis_dimension(dimension),
	m_viscosity(viscosity),
	m_dt(dt),
	m_domain(domain),
	m_num_basis_functions(dimension* dimension),
	m_prev_basis_dimension(dimension),
	m_wCoefficients(dimension* dimension, arma::fill::zeros),
	m_dwCoefficients(dimension* dimension, arma::fill::zeros),
	m_dwForceCoefficients(dimension* dimension, arma::fill::zeros),
	m_rk4Qn{ arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros) },
	m_rk4Dwt{ arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros) }
{
	Initialize();
}

// TODO: Implement the destructor

void FluidSimulation::Initialize()
{
	// TODO: Implement the Initialize function
//	CreateLookupTables();
//	CreateBasisFunctions();
}

// The structure coefficient matrices are represented here as a set of sparse matrices and are essentially a non-standard inner product between the basis functions
void FluidSimulation::CreateCkMatrices()
{
	// check if the matrices have already been created
	if (m_CkMatrices != nullptr)
	{
		DestroyCkMatrices();
	}

	// Allocate memory for the matrices
	m_CkMatrices = new arma::sp_mat[m_num_basis_functions];

	for (int k = 0; k < m_num_basis_functions; k++)
	{
		m_CkMatrices[k].set_size(m_num_basis_functions, m_num_basis_functions);
	}
}

// Create the lookup tables for the basis functions
void FluidSimulation::CreateLookupTables()
{
	if (m_basis_lookup_table != nullptr)
	{
		DestroyLookupTables();
	}

	// TODO: make sure that m_num_basis_functions and m_num_basis_functions_sqroot are set before calling this function
	// store the lookup table as a list of integer grid points
	m_basis_lookup_table = new glm::ivec2[m_num_basis_functions];

	for (int k = 0; k < m_num_basis_functions; k++)
	{
		m_basis_lookup_table[k].x = -1;
		m_basis_lookup_table[k].y = -1;
	}

	// store the reverse lookup table as a list of pointers... this is a bit of a hack but it should work TODO: explain better here
	m_basis_rlookup_table = new int* [m_basis_dimension + 1];
	for (int i = 0; i <= m_basis_dimension; i++)
	{
		m_basis_rlookup_table[i] = new int[m_basis_dimension + 1];
	}

	for (int k1 = 0; k1 <= m_basis_dimension; k1++)
	{
		for (int k2 = 0; k2 <= m_basis_dimension; k2++)
		{
			m_basis_rlookup_table[k1][k2] = -1;
		}
	}
}

// Create the scalar eigenvalues
void FluidSimulation::CreateEigenvalues()
{
	if (m_eigs != nullptr)
	{
		DestroyEigenvalues();
	}

	m_eigs = new double[m_num_basis_functions];
	m_eigs_inv = new double[m_num_basis_functions];
	m_eigs_inv_root = new double[m_num_basis_functions];

	for (int k = 0; k < m_num_basis_functions; k++)
	{
		m_eigs[k] = 0.0;
		m_eigs_inv[k] = 0.0;
		m_eigs_inv_root[k] = 0.0;
	}
}

// Fill the lookup and reverse lookup tables
void FluidSimulation::FillLookupTables()
{
	CreateLookupTables();

	int idx = 0;
	for (int k1 = 0; k1 <= m_basis_dimension; k1++)
	{
		for (int k2 = 0; k2 <= m_basis_dimension; k2++)
		{
			if (k1 > m_basis_dimension || k1 < 1 || k2 > m_basis_dimension || k2 < 1)
			{
				continue;
			}

			m_basis_lookup_table[idx].x = k1;
			m_basis_lookup_table[idx].y = k2;

			m_basis_rlookup_table[k1][k2] = idx;
			idx++;
		}
	}
}

void FluidSimulation::PrecomputeDynamics()
{
	CreateEigenvalues();

	for (int k = 0; k < m_num_basis_functions; k++)
	{
		glm::ivec2 idxK;
		idxK.x = BasisLookup(k, 0);
		idxK.y = BasisLookup(k, 1);

		double eigenvalue = (double)glm::dot(idxK, idxK);

		m_eigs[k] =	eigenvalue;
		m_eigs_inv[k] = 1.0 / eigenvalue;
		m_eigs_inv_root[k] = sqrt(m_eigs_inv[k]);
	}

	for (int d1 = 0; d1 < m_num_basis_functions; d1++)
	{
		glm::ivec2 idxA;
		idxA.x = BasisLookup(d1, 0);
		idxA.y = BasisLookup(d1, 1);

		double lambda_a = -(double)(glm::dot(idxA, idxA));

		for (int d2 = 0; d2 < m_num_basis_functions; d2++)
		{
			glm::ivec2 idxB;
			idxB.x = BasisLookup(d2, 0);
			idxB.y = BasisLookup(d2, 1);

			double lambda_b = -(double)(glm::dot(idxB, idxB));
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
				int idx = BasisRLookup(antipairs[c]);

				if (idx != -1)
				{
					glm::ivec2 idxC(c, 0);
					double coef = -CalculateCoefficient(idxA, idxB, idxC) * inv_lambda_b;
					double coef2 = -(lambda_b / lambda_a) * coef;

					if (coef != 0.0)
					{
#ifdef _WIN32
						m_CkMatrices[idx](k1, k2) = coef;
#else
						m_CkMatrices[idx].insert(k1, k2) = coef;
#endif
					}
					if (coef2 != 0.0)
					{
#ifdef _WIN32
						m_CkMatrices[idx](k2, k1) = coef2;
#else
						m_CkMatrices[idx].insert(k2, k1) = coef2;
#endif
					}					
				}
			}
		}
	}
}

//// Create the basis functions field
//void FluidSimulation::CreateBasisField()
//{
//	if (m_velocity_basis != nullptr)
//	{
//		DestroyBasisField();
//	}
//
//	// TODO: make sure that m_num_basis_functions and m_num_basis_functions_sqroot are set before calling this function and handle zero modes
//
//	// Allocate memory for the basis functions
//	m_velocity_basis = new glm::dvec2** [m_num_basis_functions];
//
//	for (int k = 0; k < m_num_basis_functions; k++)
//	{
//		m_velocity_basis[k] = new glm::dvec2 * [m_velocity_grid_resolution_y + 1];
//
//		for (int row = 0; row <= m_velocity_grid_resolution_y; row++)
//		{
//			m_velocity_basis[k][row] = new glm::dvec2[m_velocity_grid_resolution_x + 1];
//			for (int col = 0; col <= m_velocity_grid_resolution_x; col++) 
//			{
//				m_velocity_basis[k][row][col] = glm::dvec2(0.0, 0.0);
//			}
//		}
//	}		
//}

// Destroy the structure coefficient matrices
void FluidSimulation::DestroyCkMatrices()
{
	if (m_CkMatrices != nullptr)
	{
		delete[] m_CkMatrices;
		m_CkMatrices = nullptr;
	}
}

// Destroy the lookup tables
void FluidSimulation::DestroyLookupTables()
{
	if(m_basis_lookup_table != nullptr)
	{
		delete[] m_basis_lookup_table;
		m_basis_lookup_table = nullptr;
	}

	if (m_basis_rlookup_table != nullptr)
	{
		for (int i = 0; i <= m_basis_dimension; i++)
		{
			delete[] m_basis_rlookup_table[i];
			m_basis_rlookup_table[i] = nullptr;
		}
		delete[] m_basis_rlookup_table;
		m_basis_rlookup_table = nullptr;
	}
}

// Destroy the scalar eigenvalues
void FluidSimulation::DestroyEigenvalues()
{
	if (m_eigs != nullptr)
	{
		delete[] m_eigs;
		m_eigs = nullptr;
	}

	if (m_eigs_inv != nullptr)
	{
		delete[] m_eigs_inv;
		m_eigs_inv = nullptr;
	}

	if (m_eigs_inv_root != nullptr)
	{
		delete[] m_eigs_inv_root;
		m_eigs_inv_root = nullptr;
	}
}

//// Destroy the basis functions field
//void FluidSimulation::DestroyBasisField()
//{
//	if (m_velocity_basis != nullptr)
//	{
//		for (int k = 0; k < m_num_basis_functions; k++)
//		{
//			if (m_velocity_basis[k] != nullptr)
//			{
//				for (int row = 0; row <= m_velocity_grid_resolution_y; row++) 
//				{
//					if (m_velocity_basis[k][row] != nullptr)
//					{
//						delete[] m_velocity_basis[k][row];
//						m_velocity_basis[k][row] = nullptr;
//					}
//				}
//
//				delete[] m_velocity_basis[k];
//				m_velocity_basis[k] = nullptr;
//			}			
//		}
//
//		delete[] m_velocity_basis;
//		m_velocity_basis = nullptr;		
//	}
//}

// Set the basis dimension
void FluidSimulation::SetBasisDimension(int dimension)
{
	try {
		if (dimension <= 0)
		{
			throw std::invalid_argument("The basis dimension must be greater than 0.");
		}

		// Debugging //
#ifndef NDEBUG
		std::cout << "Setting the basis dimension to " << dimension << std::endl;
#endif // !NDEBUG

		// Do some teardown
		DestroyLookupTables();
		DestroyCkMatrices();
//		DestroyBasisField();
		DestroyEigenvalues();

		m_prev_basis_dimension = m_basis_dimension;

		m_basis_dimension = dimension;

		InitFields();
		FillLookupTables();
		PrecomputeBasisField();
		ExpandBasis();


	}
	catch (const std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	}
}

uint32_t FluidSimulation::BasisLookup(uint32_t idx, uint32_t component)
{
	return m_basis_lookup_table[idx][component];
}

glm::ivec2 FluidSimulation::BasisLookupK(uint32_t idx)
{
	return m_basis_lookup_table[idx];
}

uint32_t FluidSimulation::BasisRLookup(const glm::ivec2& K)
{
	if (K.x > m_basis_dimension || K.x < 1 || K.y > m_basis_dimension || K.y < 1)
	{
		return -1;
	}

	return m_basis_rlookup_table[K.x][K.y];
}

double FluidSimulation::CalculateCoefficient(const glm::ivec2& a, const glm::ivec2& b, const glm::ivec2& c)
{
	switch (c.y)
	{
	case 0:
		switch (c.x)
		{
		case 0:
			return -(double)(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return (double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return -(double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return (double)(a.x * b.y - b.x * a.y) * 0.25);
			break;
		}
		break;
	case 1:
		switch (c.x)
		{
		case 0:
			return -(double)(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return -(double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return (double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return (double)(a.x * b.y - b.x * a.y) * 0.25;
			break;
		}
		break;
	case 2:
		switch (c.x)
		{
		case 0:
			return -(double)(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return -(double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return (double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return (double)(a.x * b.y - b.x * a.y) * 0.25;
			break;
		}
		break;
	case 3:
		switch (c.x)
		{
		case 0:
			return -(double)(a.x * b.y - b.x * a.y) * 0.25;
			break;
		case 1:
			return -(double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 2:
			return (double)(a.x * b.y + b.x * a.y) * 0.25;
			break;
		case 3:
			return (double)(a.x * b.y - b.x * a.y) * 0.25;
			break;
		}
		break;
	}

	return 0;
}

double FluidSimulation::CurrentEnergy()
{
	double energy = 0.0;
	for (int k = 0; k < m_num_basis_functions; k++)
	{
		double c2 = m_wCoefficients(k) * m_wCoefficients(k);
		energy += (m_eigs_inv[k] * c2);
	}

	return energy;
}

void FluidSimulation::SetEnergy(double desired_energy)
{
	double current_energy = CurrentEnergy();
	double factor = sqrt(desired_energy) / sqrt(current_energy);

	m_wCoefficients = m_wCoefficients * factor;
}

void FluidSimulation::BasisFieldRect2D(const glm::ivec2& K, uint32_t numGridCols, uint32_t numGridRows, double domainWidth, double domainHeight, double amplitude, glm::dvec2** velocityBasisElement)
{
	double scaleX = 1.0;
	double scaleY = 1.0;

	if (K.x != 0)
	{
		scaleX = -1.0 / ((double)(glm::dot(K, K)));
	}
	if (K.y != 0)
	{
		scaleY = -1.0 / ((double)(glm::dot(K, K)));
	}

	double dx = domainWidth / (double)numGridCols;
	double dy = domainHeight / (double)numGridRows;
	double k1 = (double)K.x;
	double k2 = (double)K.y;

	for (int row = 0; row <= numGridRows; row++)
	{
		double y0 = (double)row * dy;
		double y1 = y0 + dy * 0.5;
		double sinY0K2 = sin(y0 * k2);
		double cosY1K2 = cos(y1 * k2);

		for (int col = 0; col <= numGridCols; col++)
		{
			double x0 = (double)col * dx;
			double x1 = x0 + dx * 0.5;
			double sinX0K1 = sin(x0 * k1);
			double cosX1K1 = cos(x1 * k1);

			velocityBasisElement[row][col].x = -amplitude * scaleX * k2 * sinX0K1 * cosY1K2;
			velocityBasisElement[row][col].y = amplitude * scaleY * k1 * cosX1K1 * sinY0K2;
		}
	}
}

//// Set the velocity resolution
//void FluidSimulation::SetVelocityResolution(int resolution)
//{
//	try {
//		if (resolution <= 0)
//		{
//			throw std::invalid_argument("The velocity resolution must be greater than 0.");
//		}
//
//		// Debugging //
//#ifndef NDEBUG
//		std::cout << "Setting the velocity resolution to " << resolution << std::endl;
//#endif // !NDEBUG
//
//		DestroyLookupTables();
//		DestroyCkMatrices();
//		DestroyBasisField();
//		DestroyEigenvalues();
//
//		m_velocity_grid_resolution_x = resolution;
//		m_velocity_grid_resolution_y = resolution;
//
//		InitFields();
//		FillLookupTables();
//		PrecomputeBasisField();
//		PrecomputeDynamics();
//		ExpandBasis();
//
//	}
//	catch (const std::invalid_argument& e) {
//		std::cerr << e.what() << std::endl;
//	}
//
//}