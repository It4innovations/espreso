
#include "feti4i.h"

#include "apidataprovider.h"

#include <vector>

using namespace espreso;

int main(int argc, char **argv)
{
	/* -----------------------*/
	// initialize the library //
	/* -----------------------*/
	MPI_Init(&argc, &argv);
	FETI4IInit(MPI_COMM_WORLD, 1);

	/* ---------------------------------------------------------------------------------------------------------*/
	// set options (see 'src/config/ecf/linearsolver/feti.h' for all available options of the given parameters) //
	/* ---------------------------------------------------------------------------------------------------------*/
	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts);
	FETI4ISetDefaultRealOptions(ropts);

	/* SET CONFIGURATION (here we keep the default settings)
	iopts[FETI4I_SUBDOMAINS] = 4;
	iopts[FETI4I_MAX_ITERATIONS] = 1000;
	iopts[FETI4I_FETI_METHOD] = 0; // 0: TOTAL FETI; 1: HYBRID FETI
	iopts[FETI4I_PRECONDITIONER] = 3; // 0: NONE; 1: LUMPED; 2: WEIGHT; 3: DIRICHLET
	iopts[FETI4I_CGSOLVER] = 0; // 0: PCG; 1: Pipelined PCG; 2: orthogonal PCG; 3: GMRES; 4: BICGSTAB

	ropts[FETI4I_PRECISION] = 1e-8;
	*/

	{ // this test uses the provider to compute stiffness by espreso and compute results by API

	espreso::APIDataProvider provider;
	provider.prepare();

	/* ----------------------*/
	// fill stiffness matrix //
	/* ----------------------*/
	std::vector<FETI4IInt> l2g;
	provider.fillL2G(l2g);

	FETI4IMatrix matrix;
	FETI4ICreateStiffnessMatrix(
			&matrix,
			provider.nodesSize(), l2g.data(),
			0, provider.matrixType(),
			provider.DOFs(), FETI4I_ORDER_GROUPED_DOFS);
	provider.fillMatrix([&] (FETI4IInt type, FETI4IInt size, FETI4IInt *nodes, FETI4IReal *values) {
		FETI4IAddElement(matrix, type, size, nodes, values);
	});

	/* ------------------------------------------------------------------------*/
	// create instance of the system (it also call FETI solver pre-processing) //
	/* ------------------------------------------------------------------------*/
	FETI4IInstance            instance;
	std::vector<FETI4IInt>    dirichlet_indices;
	std::vector<FETI4IReal>   dirichlet_values;
	std::vector<FETI4IMPIInt> neighbors;

	provider.fillDirichlet(dirichlet_indices, dirichlet_values);
	provider.fillNeighbors(neighbors);

	FETI4ICreateInstance(
			&instance, matrix,
			neighbors.size(), neighbors.data(), // neighbors clusters
			dirichlet_indices.size(), dirichlet_indices.data(), dirichlet_values.data(), // Dirichlet boundary condition
			iopts, ropts); // FETI4I options

	/* ---------------------------------*/
	// compute RHS and solve the system //
	/* ---------------------------------*/
	std::vector<double> rhs, solution;

	provider.fillRHS(rhs);
	solution.resize(rhs.size());

	FETI4ISolve(instance, rhs.data(), solution.data());

	/* ----------------------------------*/
	// use provider to store the results //
	/* ----------------------------------*/
	provider.storeSolution(solution);

	FETI4IDestroy(matrix);
	FETI4IDestroy(instance);

	} // destroy provider

	FETI4IFinalize();
	MPI_Finalize();
	return 0;
}


