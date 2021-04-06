
#include "include/feti4i.h"

#include "apitest.h"

int main(int argc, char** argv)
{
	// Always initialize MPI before call ESPRESO!!
	MPI_Init(&argc, &argv);

	// API Test use ESPRESO as data source for API
	//
	// Calling of provider should be replaced by appropriate methods
	espreso::APITestESPRESODataProvider provider(&argc, &argv);

	// Solving the problem by FETI4I works in 4 steps:
	//  1. Create the stiffness matrix
	//  2. Configure the ESPRESO library
	//  3. Create an instance of a problem
	//  4. Solve the instance


	// Step 1: Create the stiffness matrix.
	//       : Fill the empty matrix by element data
	FETI4IMatrix K;
	FETI4IInt    matrixType = provider.matrixType();
	FETI4IInt    indexBase = 0;

	FETI4ICreateStiffnessMatrix(&K, matrixType, indexBase); // create matrix
	for (size_t e = 0; e < provider.elements(); e++) {
		provider.addElementMatrix(K, e); // add element data
	}


	// Step 2: Configure the ESPRESO library
	//       : Set all options to default values
	//       : Change required parameters
	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts); // set default integer options
	FETI4ISetDefaultRealOptions(ropts); // set default real options

	/* here no options are changed
	iopts[FETI4IIntegerOptions::FETI4I_FETI_METHOD] = 1; // set HYBRID FETI
	iopts[FETI4IIntegerOptions::FETI4I_PRECONDITIONER] = 3; // set Dirichlet preconditioner
	*/


	// Step 3: Create an instance of a problem
	//       : Compute RHS
	//       : Fill L2G, Dirichlet, and neighbors
	FETI4IInstance            instance;
	std::vector<FETI4IReal>   rhs;
	std::vector<FETI4IInt>    dirichlet_indices;
	std::vector<FETI4IReal>   dirichlet_values;
	std::vector<FETI4IInt>    l2g;
	std::vector<FETI4IMPIInt> neighbors;

	provider.computeRHS(rhs);
	provider.fillL2G(l2g);
	provider.fillDirichlet(dirichlet_indices, dirichlet_values);
	provider.fillNeighbours(neighbors);

	FETI4ICreateInstance(
			&instance,
			K, // Stiffness matrix
			rhs.size(), rhs.data(), // RHS
			l2g.data(), // local 2 global indices mapping
			neighbors.size(), neighbors.data(), // neighbors clusters
			dirichlet_indices.size(), dirichlet_indices.data(), dirichlet_values.data(), // Dirichlet boundary condition
			iopts, ropts); // FETI4I options


	// Step 4: Solve the instance
	//       : Process the solution
	std::vector<FETI4IReal> solution(rhs.size());
	FETI4ISolve(instance, solution.size(), solution.data());


	// Finish: Destroy all data
	FETI4IDestroy(K);
	FETI4IDestroy(instance);

	MPI_Finalize();
}

