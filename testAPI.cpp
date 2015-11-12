
#include "libespreso/espreso.h"
#include <sstream>
#include <string>
#include <iostream>

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	int MPIrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	ESPRESOInit(MPI_COMM_WORLD);

	ESPRESOMat K;
	//ESPRESOCreateMatrixElemental(n, nelt, eltptr, eltvar, values, &K);

	std::stringstream ss;
	ss << argv[1] << MPIrank;

	std::stringstream KeInfo;

	KeInfo << ss.str() << "KeInfo.txt";

	std::cout << ss.str() << "\n";



	ESPRESOFinalize();

	MPI_Finalize();
}



