
#include "../configuration/environment.h"

#include "../basis/utilities/utils.h"

namespace espreso {

Environment *environment;

Environment::Environment(): executable("espreso")
{
	int initialized;
	MPI_Initialized(&initialized);
	if (initialized) {
		MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
		MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	}

	executable = Esutils::getEnv<std::string>("_");

	MKL_NUM_THREADS    = Esutils::getEnv<size_t>("MKL_NUM_THREADS");
	OMP_NUM_THREADS    = Esutils::getEnv<size_t>("OMP_NUM_THREADS");
	SOLVER_NUM_THREADS = Esutils::getEnv<size_t>("SOLVER_NUM_THREADS");
	PAR_NUM_THREADS    = Esutils::getEnv<size_t>("PAR_NUM_THREADS");
	CILK_NWORKERS      = Esutils::getEnv<size_t>("CILK_NWORKERS");

	environment = this;
}

}



