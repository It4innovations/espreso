
#include "factory/factory.h"

using namespace espreso;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &config::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::MPIsize);

	Options options(&argc, &argv);

	//std::cout << "USE_SCHUR_COMPLEMENT - " << config::solver::USE_SCHUR_COMPLEMENT << std::endl;
	//std::cout << "SCHUR_COMPLEMENT_PREC - " << config::solver::SCHUR_COMPLEMENT_PREC << std::endl;
	//std::cout << "COMBINE_SC_AND_SPDS - " << config::solver::COMBINE_SC_AND_SPDS << std::endl;

	Factory factory(options);
	factory.solve(1);
	factory.store("result");

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


