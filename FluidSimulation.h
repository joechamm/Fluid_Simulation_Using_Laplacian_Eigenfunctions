#pragma once

#define _USE_MATH_DEFINES
#include <glm/glm.hpp>
#include <armadillo>
#include <cmath>

/*
	The fluid simulation class will run our simulation and store the state of the simulation. The main objects of our simulation are the eigenfunctions of the 
	vector Laplacian operator. We're assuming a 2D rectangular domain [0,pi]x[0,pi] here, so the eigenfunctions are just sin's and cos's. We'll index the basis 
	functions by k = (k1, k2) where k1 and k2 are integers. This gives us a closed form expression for the k'th basis function:
		phi_k = (1/(k1^2 + k2^2)) * <k2 * sin(k1 * x) * cos(k2 * y), -k1 * cos(k1 * x) * sin(k2 * y)>
	The vector part is <k2 * sin(k1 * x) * cos(k2 * y), -k1 * cos(k1 * x) * sin(k2 * y)> and the associated eigenvalue is - (k1^2 + k2^2). 
	Since we have a closed form expression for the basis functions, we don't need to store them explicitly. If we need to evaluate the velocity, we will sample
	it. We just need to store the coefficients of the basis functions, which uniquely determine the state of the simulation. We'll precompute the eigenvalues,
	multiplicative inverses of the eigenvalues, and the square root of the multiplicative inverses of the eigenvalues for efficiency. We'll also precompute the
	structure coefficient matrices, which allow us a non-standard inner product between pairs of basis functions. This is needed for the advection step of the
	simulation which abstractly can be thought of here as the Lie derivative, or Jacobi-Lie bracket of vector fields -[u,w], where u is the velocity field and w
	is the vorticity field. Since the Jacobi-Lie bracket and vector cross product are anti-symmetric operators, the structure coefficient matrices have the property 
	that (1/lamba_i)C_k[i,j] = (-1/lambda_j)C_k[j,i]. This allows us to compute the advection step efficiently. We'll also precompute the basis function lookup
	tables and reverse lookup tables for efficiency. We also need to store the previous state of the simulation for the advection step. We'll use a 4th order
	Runge-Kutta integration scheme for the advection step. We'll also store the external force projection coefficients for the external force projection step of the
	simulation. We'll also store the zero modes of the basis functions, which we can use to see the effect of different basis functions. 
*/

class FluidSimulation
{
	// Basis function lookup table
	glm::ivec2* m_basis_lookup_table;
	// Reverse lookup table
	int32_t** m_basis_rlookup_table;
	// Scalar Eigenvalues
	double* m_eigs;
	// Multiplicative inverses of eigenvalues
	double* m_eigs_inv;
	// Square root of multiplicative inverses of eigenvalues
	double* m_eigs_inv_root;

	// Structure Coefficient Matrices
	arma::sp_mat* m_CkMatrices;
	// Basis coefficients
	arma::vec m_wCoefficients;
	// Basis coefficients delta for timestep integration
	arma::vec m_dwCoefficients;
	// External force projection coefficients
	arma::vec m_dwForceCoefficients;
	// Basis coefficients for 4th order RungaKutta integration
	arma::vec m_rk4Qn[4];
	arma::vec m_rk4Dwt[4];

	/// Simulation parameters ///
	// viscosity
	double m_viscosity;
	// timestep
	double m_dt;

	// simulation domain, probably just [0,pi], [0,pi] for now
	// store in a 4D vector (x_min, y_min, x_max, y_max), so for example [0,pi]x[0,pi] would be (0,0,pi,pi)
	glm::dvec4 m_domain;  

	// since we are using a square grid, we can use the square root of the number of basis functions, which is the number of basis functions in one dimension (though it doesn't quite work the same way here since our field represents the curl)
	// we'll call this the dimension of the simulation and it will be the square root of the number of basis functions
	uint32_t m_basis_dimension; // just call it dimension, then the number of basis functions is dimension^2
	// number of basis functions
	uint32_t m_num_basis_functions;

public:
	FluidSimulation(uint32_t dimension = 6, double viscosity = 0.01, double dt = 0.1, const glm::dvec4& domain = glm::dvec4(0.0, 0.0, M_PI, M_PI));
	~FluidSimulation();

	/// Accessors
	double GetViscosity() const { return m_viscosity; }
	void SetViscosity(double viscosity) { m_viscosity = viscosity; }
	double GetTimestep() const { return m_dt; }
	void SetTimestep(double dt) { m_dt = dt; }
	glm::dvec4 GetDomain() const { return m_domain; }
	void SetDomain(const glm::dvec4& domain) { m_domain = domain; }
	uint32_t GetBasisDimension() const { return m_basis_dimension; }
	uint32_t GetNumBasisFunctions() const { return m_num_basis_functions; }
	// calculate the energy of the simulation
	double CurrentEnergy() const;
	// set the energy of the simulation
	void SetEnergy(double desired_energy);

	// Basis lookup table accessors
	int32_t BasisLookup(uint32_t idx, uint32_t component) const;
	glm::ivec2 BasisLookupK(uint32_t idx) const;
	int32_t BasisRLookup(const glm::ivec2& K) const;

	// Calculate the coefficient
	double CalculateCoefficient(const glm::ivec2& a, const glm::ivec2& b, const glm::ivec2& c) const;

	// project the list of external forces onto the basis functions where the first two components are the x and y position components and the last two are the x and y components of the force vector
	void ProjectForces(const std::vector<glm::dvec4>& force_list, arma::vec& delW) const;
	// calculate the velocity basis element for a given basis function
	void BasisFieldRect2D(const glm::ivec2& K, uint32_t numGridCols, uint32_t numGridRows, double amplitude, glm::dvec2** velocityBasisElement) const;

	// step the simulation forward in time
	void TimeStep();

	// Initialize the simulation
	void Initialize();

	// Fill Lookup Table
	void FillLookupTables();

	void PrecomputeDynamics();
	
	/// Teardown ///
	// Destroy the structure coefficient matrices
	void DestroyCkMatrices();
	// Destroy the lookup tables
	void DestroyLookupTables();
	// Destroy the scalar eigenvalues
	void DestroyEigenvalues();
};

