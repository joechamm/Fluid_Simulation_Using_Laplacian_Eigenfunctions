#include "../../include/FluidSimulation.h"

// TODO: set m_velocity_grid_resolution_y and m_velocity_grid_resolution_x 
FluidSimulation::FluidSimulation(uint32_t dimension, double viscosity, double dt, const glm::dvec4& domain) :
	m_basis_lookup_table(nullptr),
	m_basis_rlookup_table(nullptr),
	m_eigs(nullptr),
	m_eigs_inv(nullptr),
	m_eigs_inv_root(nullptr),
	m_CkMatrices(nullptr),
	m_basis_dimension(dimension),
	m_viscosity(viscosity),
	m_dt(dt),
	m_domain(domain),
	m_num_basis_functions(dimension* dimension),
	m_wCoefficients(dimension* dimension, arma::fill::zeros),
	m_dwCoefficients(dimension* dimension, arma::fill::zeros),
	m_dwForceCoefficients(dimension* dimension, arma::fill::zeros),
	m_rk4Qn{ arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros) },
	m_rk4Dwt{ arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros), arma::vec(dimension* dimension, arma::fill::zeros) }
{
	// create the lookup tables
	m_basis_lookup_table = new glm::ivec2[m_num_basis_functions];

	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		m_basis_lookup_table[k].x = -1;
		m_basis_lookup_table[k].y = -1;
	}

	m_basis_rlookup_table = new int32_t * [m_basis_dimension + 1];
	for (uint32_t i = 0; i <= m_basis_dimension; i++)
	{
		m_basis_rlookup_table[i] = new int32_t[m_basis_dimension + 1];
	}

	for(uint32_t k1 = 0; k1 <= m_basis_dimension; k1++)
	{
		for (uint32_t k2 = 0; k2 <= m_basis_dimension; k2++)
		{
			m_basis_rlookup_table[k1][k2] = -1;
		}
	}

	// create the eigenvalues
	m_eigs = new double[m_num_basis_functions];
	m_eigs_inv = new double[m_num_basis_functions];
	m_eigs_inv_root = new double[m_num_basis_functions];

	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		m_eigs[k] = 0.0;
		m_eigs_inv[k] = 0.0;
		m_eigs_inv_root[k] = 0.0;
	}

	// create the structure coefficient matrices
	m_CkMatrices = new arma::sp_mat[m_num_basis_functions];

	for (uint32_t k; k < m_num_basis_functions; k++)
	{
		m_CkMatrices[k].set_size(m_num_basis_functions, m_num_basis_functions);
	}
}

// TODO: Implement the destructor
FluidSimulation::~FluidSimulation()
{
	DestroyLookupTables();
	DestroyEigenvalues();
	DestroyCkMatrices();
}

void FluidSimulation::Initialize()
{
	// TODO: Implement the Initialize function
//	CreateLookupTables();
//	CreateBasisFunctions();
	m_wCoefficients(0) = 1.0;
	FillLookupTables();
	// PrecomputeBasisFields();
	PrecomputeDynamics();
}

// Fill the lookup and reverse lookup tables
void FluidSimulation::FillLookupTables()
{
	int32_t idx = 0;
	for (uint32_t k1 = 0; k1 <= m_basis_dimension; k1++)
	{
		for (uint32_t k2 = 0; k2 <= m_basis_dimension; k2++)
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
	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		glm::ivec2 idxK;
		idxK.x = BasisLookup(k, 0);
		idxK.y = BasisLookup(k, 1);

		double eigenvalue = (double)glm::dot(idxK, idxK);

		m_eigs[k] =	eigenvalue;
		m_eigs_inv[k] = 1.0 / eigenvalue;
		m_eigs_inv_root[k] = sqrt(m_eigs_inv[k]);
	}

	for (uint32_t d1 = 0; d1 < m_num_basis_functions; d1++)
	{
		glm::ivec2 idxA;
		idxA.x = BasisLookup(d1, 0);
		idxA.y = BasisLookup(d1, 1);

		double lambda_a = -(double)(glm::dot(idxA, idxA));

		for (uint32_t d2 = 0; d2 < m_num_basis_functions; d2++)
		{
			glm::ivec2 idxB;
			idxB.x = BasisLookup(d2, 0);
			idxB.y = BasisLookup(d2, 1);

			double lambda_b = -(double)(glm::dot(idxB, idxB));
			double inv_lambda_b = 1.0 / lambda_b;

			int32_t k1 = BasisRLookup(idxA);
			int32_t k2 = BasisRLookup(idxB);

			glm::ivec2 antipairs[4];

			antipairs[0].x = idxA.x - idxB.x;
			antipairs[1].x = idxA.x - idxB.x;
			antipairs[2].x = idxA.x + idxB.x;
			antipairs[3].x = idxA.x + idxB.x;

			antipairs[0].y = idxA.y - idxB.y;
			antipairs[1].y = idxA.y + idxB.y;
			antipairs[2].y = idxA.y - idxB.y;
			antipairs[3].y = idxA.y + idxB.y;

			for (uint32_t c = 0; c < 4; c++)
			{
				int32_t idx = BasisRLookup(antipairs[c]);

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

int32_t FluidSimulation::BasisLookup(uint32_t idx, uint32_t component) const
{
	return m_basis_lookup_table[idx][component];
}

glm::ivec2 FluidSimulation::BasisLookupK(uint32_t idx) const
{
	return m_basis_lookup_table[idx];
}

int32_t FluidSimulation::BasisRLookup(const glm::ivec2& K) const
{
	if (K.x > m_basis_dimension || K.x < 1 || K.y > m_basis_dimension || K.y < 1)
	{
		return -1;
	}

	return m_basis_rlookup_table[K.x][K.y];
}

double FluidSimulation::CalculateCoefficient(const glm::ivec2& a, const glm::ivec2& b, const glm::ivec2& c) const
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
			return (double)(a.x * b.y - b.x * a.y) * 0.25;
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

double FluidSimulation::CurrentEnergy() const
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

void FluidSimulation::BasisFieldRect2D(const glm::ivec2& K, uint32_t numGridCols, uint32_t numGridRows, double amplitude, glm::dvec2** velocityBasisElement) const
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

	double domainWidth = m_domain.z - m_domain.x;
	double domainHeight = m_domain.w - m_domain.y;

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

void FluidSimulation::TimeStep()
{
	// store the current energy for renomalization later
	double previous_energy = CurrentEnergy();

	// 4th order-Runge-Kutta integration
	m_rk4Qn[0] = m_wCoefficients;

	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		m_rk4Dwt[0](k) = arma::as_scalar(arma::trans(m_rk4Qn[0]) * m_CkMatrices[k] * m_rk4Qn[0]);
		m_rk4Qn[1](k) = m_rk4Qn[0](k) + m_rk4Dwt[0](k) * m_dt * 0.5;
	}

	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		m_rk4Dwt[1](k) = arma::as_scalar(arma::trans(m_rk4Qn[1]) * m_CkMatrices[k] * m_rk4Qn[1]);
		m_rk4Qn[2](k) = m_rk4Qn[0](k) + m_rk4Dwt[1](k) * m_dt * 0.5;
	}

	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		m_rk4Dwt[2](k) = arma::as_scalar(arma::trans(m_rk4Qn[2]) * m_CkMatrices[k] * m_rk4Qn[2]);
		m_rk4Qn[3](k) = m_rk4Qn[0](k) + m_rk4Dwt[2](k) * m_dt;
	}

	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		m_rk4Dwt[3](k) = arma::as_scalar(arma::trans(m_rk4Qn[3]) * m_CkMatrices[k] * m_rk4Qn[3]);
		m_dwCoefficients(k) = (m_rk4Dwt[0](k) + m_rk4Dwt[1](k) * 2.0 + m_rk4Dwt[2](k) * 2.0 + m_rk4Dwt[3](k)) / 6.0;
	}

	m_wCoefficients += m_dwCoefficients * m_dt;

	if (previous_energy > 1e-5)
	{
		SetEnergy(previous_energy);
	}

	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		double eigenvalue = - m_eigs[k];
		m_wCoefficients(k) = m_wCoefficients(k) * exp(eigenvalue * m_dt * m_viscosity) + m_dwForceCoefficients(k);
		m_dwForceCoefficients(k) = 0.0;
	}
}

void FluidSimulation::ProjectForces(const std::vector<glm::dvec4>& force_list, arma::vec& delW) const
{
	delW.set_size(m_num_basis_functions);
	delW.zeros();
	for (uint32_t k = 0; k < m_num_basis_functions; k++)
	{
		int32_t k1 = BasisLookup(k, 0);
		int32_t k2 = BasisLookup(k, 1);

		double xfactor = 1.0;
		double yfactor = 1.0;
		if (k1 != 0)
		{
			xfactor = -1.0 / ((double)(k1 * k1 + k2 * k2));
		}
		if (k2 != 0)
		{
			yfactor = -1.0 / ((double)(k1 * k1 + k2 * k2));
		}

		std::vector<glm::dvec4>::const_iterator it;
		for (it = force_list.begin(); it != force_list.end(); ++it)
		{
			double x = it->x;
			double y = it->y;
			double fx = it->z;
			double fy = it->w;

			if (x >= 1.00001 || x <= -0.00001 || y >= 1.00001 || y <= -0.00001)
			{
				continue;
			}

			x *= (m_domain.z - m_domain.x);
			y *= (m_domain.w - m_domain.y);

			double vx = -((double)k2) * xfactor * sin(((double)k1) * x) * cos(((double)k2) * y) * m_dt;
			double vy = ((double)k1) * yfactor * cos(((double)k1) * x) * sin(((double)k2) * y) * m_dt;

			delW(k) += (vx * fx + vy * fy);
		}
	}
}
