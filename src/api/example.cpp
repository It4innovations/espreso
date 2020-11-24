
#include "feti4i.h"

#include <vector>

std::vector<FETI4IReal> stiffness();

int main(int argc, char **argv)
{
	/* -----------------------*/
	// initialize the library //
	/* -----------------------*/
	MPI_Init(&argc, &argv);
	FETI4IInit(MPI_COMM_WORLD, 0);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (size != 2) {
		printf("RUN THIS EXAMPLE WITH 2 MPI PROCESSES ONLY!\n");
		return 1;
	}

	/* ---------------------------------------------------------------------------------------------------------*/
	// set options (see 'src/config/ecf/linearsolver/feti.h' for all available options of the given parameters) //
	/* ---------------------------------------------------------------------------------------------------------*/
	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts);
	FETI4ISetDefaultRealOptions(ropts);

	iopts[FETI4I_SUBDOMAINS] = 1;
	/* SET CONFIGURATION (here we keep the default settings)
	iopts[FETI4I_MAX_ITERATIONS] = 1000;
	iopts[FETI4I_FETI_METHOD] = 0; // 0: TOTAL FETI; 1: HYBRID FETI
	iopts[FETI4I_PRECONDITIONER] = 3; // 0: NONE; 1: LUMPED; 2: WEIGHT; 3: DIRICHLET
	iopts[FETI4I_CGSOLVER] = 0; // 0: PCG; 1: Pipelined PCG; 2: orthogonal PCG; 3: GMRES; 4: BICGSTAB

	ropts[FETI4I_PRECISION] = 1e-8;
	*/

	{ // this test uses the provider to compute stiffness by espreso and compute results by API

	/* ----------------------*/
	// fill stiffness matrix //
	/* ----------------------*/
	FETI4IInt nnodes = 8, dofs = 3, mtype = FETI4I_REAL_SYMMETRIC_POSITIVE_DEFINITE;
	std::vector<FETI4IInt> l2g, nodes =  { 0, 1, 3, 2, 4, 5, 7, 6 };
	std::vector<FETI4IReal> values = stiffness();
	if (rank) {
		l2g = { 4, 5, 6, 7, 8, 9, 10, 11 };
	} else {
		l2g = { 0, 1, 2, 3, 4, 5, 6, 7 };
	}

	FETI4IMatrix matrix;
	FETI4ICreateStiffnessMatrix( &matrix, nnodes, l2g.data(), 0, mtype, dofs, FETI4I_ORDER_GROUPED_DOFS);
	FETI4IAddElement(matrix, FETI4I_ETYPE_VOLUME, 8, nodes.data(), values.data());

	/* ------------------------------------------------------------------------*/
	// create instance of the system (it also call FETI solver pre-processing) //
	/* ------------------------------------------------------------------------*/
	FETI4IInstance            instance;
	std::vector<FETI4IInt>    dirichlet_indices;
	std::vector<FETI4IReal>   dirichlet_values;
	std::vector<FETI4IMPIInt> neighbors = { (rank + 1) % size };
	if (rank == 0) {
		// index = node * dofs + dof
		dirichlet_indices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
		dirichlet_values = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	}

	FETI4ICreateInstance(
			&instance, matrix,
			neighbors.size(), neighbors.data(), // neighbors clusters
			dirichlet_indices.size(), dirichlet_indices.data(), dirichlet_values.data(), // Dirichlet boundary condition
			iopts, ropts); // FETI4I options

	/* ---------------------------------*/
	// compute RHS and solve the system //
	/* ---------------------------------*/
	std::vector<double> rhs = { 0, 0, 4.81136e+09, 0, 0, 4.81136e+09, 0, 0, 4.81136e+09, 0, 0, 4.81136e+09, 0, 0, 4.81136e+09, 0, 0, 4.81136e+09, 0, 0, 4.81136e+09, 0, 0, 4.81136e+09 };
	std::vector<double> solution(rhs.size());
	FETI4ISolve(instance, rhs.data(), solution.data());

	FETI4IDestroy(matrix);
	FETI4IDestroy(instance);

	} // destroy provider

	FETI4IFinalize();
	MPI_Finalize();
	return 0;
}

std::vector<FETI4IReal> stiffness()
{
	return {
			3.8141e+10, -4.48718e+09, -5.60897e+09, 1.23397e+10, -7.85256e+09, -1.57051e+10, -9.53526e+09, -7.29167e+09, 8.41346e+09, 1.68269e+09, -8.41346e+09, -1.68269e+09, 4.20673e+09, 8.41346e+08, -4.20673e+09, -8.41346e+08, 1.68269e+10, 3.36538e+09, 1.68269e+09, 8.41346e+09, -3.36538e+09, -1.68269e+10, -8.41346e+09, -1.68269e+09,
			-4.48718e+09, 3.8141e+10, 1.23397e+10, -5.60897e+09, -1.57051e+10, -7.85256e+09, -7.29167e+09, -9.53526e+09, -1.68269e+09, -8.41346e+09, 1.68269e+09, 8.41346e+09, -8.41346e+08, -4.20673e+09, 8.41346e+08, 4.20673e+09, -3.36538e+09, -1.68269e+10, -8.41346e+09, -1.68269e+09, 1.68269e+10, 3.36538e+09, 1.68269e+09, 8.41346e+09,
			-5.60897e+09, 1.23397e+10, 3.8141e+10, -4.48718e+09, -9.53526e+09, -7.29167e+09, -7.85256e+09, -1.57051e+10, -8.41346e+09, -1.68269e+09, 8.41346e+09, 1.68269e+09, -4.20673e+09, -8.41346e+08, 4.20673e+09, 8.41346e+08, -1.68269e+09, -8.41346e+09, -1.68269e+10, -3.36538e+09, 8.41346e+09, 1.68269e+09, 3.36538e+09, 1.68269e+10,
			1.23397e+10, -5.60897e+09, -4.48718e+09, 3.8141e+10, -7.29167e+09, -9.53526e+09, -1.57051e+10, -7.85256e+09, 1.68269e+09, 8.41346e+09, -1.68269e+09, -8.41346e+09, 8.41346e+08, 4.20673e+09, -8.41346e+08, -4.20673e+09, 8.41346e+09, 1.68269e+09, 3.36538e+09, 1.68269e+10, -1.68269e+09, -8.41346e+09, -1.68269e+10, -3.36538e+09,
			-7.85256e+09, -1.57051e+10, -9.53526e+09, -7.29167e+09, 3.8141e+10, -4.48718e+09, -5.60897e+09, 1.23397e+10, 4.20673e+09, 8.41346e+08, -4.20673e+09, -8.41346e+08, 8.41346e+09, 1.68269e+09, -8.41346e+09, -1.68269e+09, 3.36538e+09, 1.68269e+10, 8.41346e+09, 1.68269e+09, -1.68269e+10, -3.36538e+09, -1.68269e+09, -8.41346e+09,
			-1.57051e+10, -7.85256e+09, -7.29167e+09, -9.53526e+09, -4.48718e+09, 3.8141e+10, 1.23397e+10, -5.60897e+09, -8.41346e+08, -4.20673e+09, 8.41346e+08, 4.20673e+09, -1.68269e+09, -8.41346e+09, 1.68269e+09, 8.41346e+09, -1.68269e+10, -3.36538e+09, -1.68269e+09, -8.41346e+09, 3.36538e+09, 1.68269e+10, 8.41346e+09, 1.68269e+09,
			-9.53526e+09, -7.29167e+09, -7.85256e+09, -1.57051e+10, -5.60897e+09, 1.23397e+10, 3.8141e+10, -4.48718e+09, -4.20673e+09, -8.41346e+08, 4.20673e+09, 8.41346e+08, -8.41346e+09, -1.68269e+09, 8.41346e+09, 1.68269e+09, -8.41346e+09, -1.68269e+09, -3.36538e+09, -1.68269e+10, 1.68269e+09, 8.41346e+09, 1.68269e+10, 3.36538e+09,
			-7.29167e+09, -9.53526e+09, -1.57051e+10, -7.85256e+09, 1.23397e+10, -5.60897e+09, -4.48718e+09, 3.8141e+10, 8.41346e+08, 4.20673e+09, -8.41346e+08, -4.20673e+09, 1.68269e+09, 8.41346e+09, -1.68269e+09, -8.41346e+09, 1.68269e+09, 8.41346e+09, 1.68269e+10, 3.36538e+09, -8.41346e+09, -1.68269e+09, -3.36538e+09, -1.68269e+10,
			8.41346e+09, -1.68269e+09, -8.41346e+09, 1.68269e+09, 4.20673e+09, -8.41346e+08, -4.20673e+09, 8.41346e+08, 3.8141e+10, 1.23397e+10, -5.60897e+09, -4.48718e+09, -7.85256e+09, -7.29167e+09, -9.53526e+09, -1.57051e+10, 1.68269e+10, 8.41346e+09, 1.68269e+09, 3.36538e+09, -3.36538e+09, -1.68269e+09, -8.41346e+09, -1.68269e+10,
			1.68269e+09, -8.41346e+09, -1.68269e+09, 8.41346e+09, 8.41346e+08, -4.20673e+09, -8.41346e+08, 4.20673e+09, 1.23397e+10, 3.8141e+10, -4.48718e+09, -5.60897e+09, -7.29167e+09, -7.85256e+09, -1.57051e+10, -9.53526e+09, 8.41346e+09, 1.68269e+10, 3.36538e+09, 1.68269e+09, -1.68269e+09, -3.36538e+09, -1.68269e+10, -8.41346e+09,
			-8.41346e+09, 1.68269e+09, 8.41346e+09, -1.68269e+09, -4.20673e+09, 8.41346e+08, 4.20673e+09, -8.41346e+08, -5.60897e+09, -4.48718e+09, 3.8141e+10, 1.23397e+10, -9.53526e+09, -1.57051e+10, -7.85256e+09, -7.29167e+09, -1.68269e+09, -3.36538e+09, -1.68269e+10, -8.41346e+09, 8.41346e+09, 1.68269e+10, 3.36538e+09, 1.68269e+09,
			-1.68269e+09, 8.41346e+09, 1.68269e+09, -8.41346e+09, -8.41346e+08, 4.20673e+09, 8.41346e+08, -4.20673e+09, -4.48718e+09, -5.60897e+09, 1.23397e+10, 3.8141e+10, -1.57051e+10, -9.53526e+09, -7.29167e+09, -7.85256e+09, -3.36538e+09, -1.68269e+09, -8.41346e+09, -1.68269e+10, 1.68269e+10, 8.41346e+09, 1.68269e+09, 3.36538e+09,
			4.20673e+09, -8.41346e+08, -4.20673e+09, 8.41346e+08, 8.41346e+09, -1.68269e+09, -8.41346e+09, 1.68269e+09, -7.85256e+09, -7.29167e+09, -9.53526e+09, -1.57051e+10, 3.8141e+10, 1.23397e+10, -5.60897e+09, -4.48718e+09, 3.36538e+09, 1.68269e+09, 8.41346e+09, 1.68269e+10, -1.68269e+10, -8.41346e+09, -1.68269e+09, -3.36538e+09,
			8.41346e+08, -4.20673e+09, -8.41346e+08, 4.20673e+09, 1.68269e+09, -8.41346e+09, -1.68269e+09, 8.41346e+09, -7.29167e+09, -7.85256e+09, -1.57051e+10, -9.53526e+09, 1.23397e+10, 3.8141e+10, -4.48718e+09, -5.60897e+09, 1.68269e+09, 3.36538e+09, 1.68269e+10, 8.41346e+09, -8.41346e+09, -1.68269e+10, -3.36538e+09, -1.68269e+09,
			-4.20673e+09, 8.41346e+08, 4.20673e+09, -8.41346e+08, -8.41346e+09, 1.68269e+09, 8.41346e+09, -1.68269e+09, -9.53526e+09, -1.57051e+10, -7.85256e+09, -7.29167e+09, -5.60897e+09, -4.48718e+09, 3.8141e+10, 1.23397e+10, -8.41346e+09, -1.68269e+10, -3.36538e+09, -1.68269e+09, 1.68269e+09, 3.36538e+09, 1.68269e+10, 8.41346e+09,
			-8.41346e+08, 4.20673e+09, 8.41346e+08, -4.20673e+09, -1.68269e+09, 8.41346e+09, 1.68269e+09, -8.41346e+09, -1.57051e+10, -9.53526e+09, -7.29167e+09, -7.85256e+09, -4.48718e+09, -5.60897e+09, 1.23397e+10, 3.8141e+10, -1.68269e+10, -8.41346e+09, -1.68269e+09, -3.36538e+09, 3.36538e+09, 1.68269e+09, 8.41346e+09, 1.68269e+10,
			1.68269e+10, -3.36538e+09, -1.68269e+09, 8.41346e+09, 3.36538e+09, -1.68269e+10, -8.41346e+09, 1.68269e+09, 1.68269e+10, 8.41346e+09, -1.68269e+09, -3.36538e+09, 3.36538e+09, 1.68269e+09, -8.41346e+09, -1.68269e+10, 7.17949e+10, 2.91667e+10, 1.12179e+10, 2.91667e+10, -5.83333e+10, -3.25321e+10, -1.79487e+10, -3.25321e+10,
			3.36538e+09, -1.68269e+10, -8.41346e+09, 1.68269e+09, 1.68269e+10, -3.36538e+09, -1.68269e+09, 8.41346e+09, 8.41346e+09, 1.68269e+10, -3.36538e+09, -1.68269e+09, 1.68269e+09, 3.36538e+09, -1.68269e+10, -8.41346e+09, 2.91667e+10, 7.17949e+10, 2.91667e+10, 1.12179e+10, -3.25321e+10, -5.83333e+10, -3.25321e+10, -1.79487e+10,
			1.68269e+09, -8.41346e+09, -1.68269e+10, 3.36538e+09, 8.41346e+09, -1.68269e+09, -3.36538e+09, 1.68269e+10, 1.68269e+09, 3.36538e+09, -1.68269e+10, -8.41346e+09, 8.41346e+09, 1.68269e+10, -3.36538e+09, -1.68269e+09, 1.12179e+10, 2.91667e+10, 7.17949e+10, 2.91667e+10, -1.79487e+10, -3.25321e+10, -5.83333e+10, -3.25321e+10,
			8.41346e+09, -1.68269e+09, -3.36538e+09, 1.68269e+10, 1.68269e+09, -8.41346e+09, -1.68269e+10, 3.36538e+09, 3.36538e+09, 1.68269e+09, -8.41346e+09, -1.68269e+10, 1.68269e+10, 8.41346e+09, -1.68269e+09, -3.36538e+09, 2.91667e+10, 1.12179e+10, 2.91667e+10, 7.17949e+10, -3.25321e+10, -1.79487e+10, -3.25321e+10, -5.83333e+10,
			-3.36538e+09, 1.68269e+10, 8.41346e+09, -1.68269e+09, -1.68269e+10, 3.36538e+09, 1.68269e+09, -8.41346e+09, -3.36538e+09, -1.68269e+09, 8.41346e+09, 1.68269e+10, -1.68269e+10, -8.41346e+09, 1.68269e+09, 3.36538e+09, -5.83333e+10, -3.25321e+10, -1.79487e+10, -3.25321e+10, 7.17949e+10, 2.91667e+10, 1.12179e+10, 2.91667e+10,
			-1.68269e+10, 3.36538e+09, 1.68269e+09, -8.41346e+09, -3.36538e+09, 1.68269e+10, 8.41346e+09, -1.68269e+09, -1.68269e+09, -3.36538e+09, 1.68269e+10, 8.41346e+09, -8.41346e+09, -1.68269e+10, 3.36538e+09, 1.68269e+09, -3.25321e+10, -5.83333e+10, -3.25321e+10, -1.79487e+10, 2.91667e+10, 7.17949e+10, 2.91667e+10, 1.12179e+10,
			-8.41346e+09, 1.68269e+09, 3.36538e+09, -1.68269e+10, -1.68269e+09, 8.41346e+09, 1.68269e+10, -3.36538e+09, -8.41346e+09, -1.68269e+10, 3.36538e+09, 1.68269e+09, -1.68269e+09, -3.36538e+09, 1.68269e+10, 8.41346e+09, -1.79487e+10, -3.25321e+10, -5.83333e+10, -3.25321e+10, 1.12179e+10, 2.91667e+10, 7.17949e+10, 2.91667e+10,
			-1.68269e+09, 8.41346e+09, 1.68269e+10, -3.36538e+09, -8.41346e+09, 1.68269e+09, 3.36538e+09, -1.68269e+10, -1.68269e+10, -8.41346e+09, 1.68269e+09, 3.36538e+09, -3.36538e+09, -1.68269e+09, 8.41346e+09, 1.68269e+10, -3.25321e+10, -1.79487e+10, -3.25321e+10, -5.83333e+10, 2.91667e+10, 1.12179e+10, 2.91667e+10, 7.17949e+10
	};
}
