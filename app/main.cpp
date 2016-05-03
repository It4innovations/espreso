
#include "factory/factory.h"

using namespace espreso;

int main(int argc, char **argv)
{
	//ESINFO(ERROR) << "test failure";

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);

	Options options(&argc, &argv);

	Factory factory(options);
	factory.solve(1);
	factory.store("result");

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


