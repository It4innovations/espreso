
#include "mpi.h"

#include "environment.h"
#include "configuration.hpp"

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

	char *var = getenv("_");
	if (var != NULL) {
		executable = Esutils::getEnv<std::string>("_");
	} else {
		executable = "espreso";
	}

	MKL_NUM_THREADS    = ParameterHolder::create<size_t>("MKL_NUM_THREADS"   , "Number of MKL threads"          , MKL_NUM_THREADS   , Esutils::getEnv<size_t>("MKL_NUM_THREADS")   , "size_t", this);
	OMP_NUM_THREADS    = ParameterHolder::create<size_t>("OMP_NUM_THREADS"   , "Number of OMP threads"          , OMP_NUM_THREADS   , Esutils::getEnv<size_t>("OMP_NUM_THREADS")   , "size_t", this);
	SOLVER_NUM_THREADS = ParameterHolder::create<size_t>("SOLVER_NUM_THREADS", "Number of ESRESO solver threads", SOLVER_NUM_THREADS, Esutils::getEnv<size_t>("SOLVER_NUM_THREADS"), "size_t", this);
	PAR_NUM_THREADS    = ParameterHolder::create<size_t>("PAR_NUM_THREADS"   , "Number of parallel threads"     , PAR_NUM_THREADS   , Esutils::getEnv<size_t>("PAR_NUM_THREADS")   , "size_t", this);
	CILK_NWORKERS      = ParameterHolder::create<size_t>("CILK_NWORKERS"     , "Number of cilk++ threads"       , CILK_NWORKERS     , Esutils::getEnv<size_t>("CILK_NWORKERS")     , "size_t", this);

	log_dir = ParameterHolder::create<std::string>("log_dir", "Log directory", log_dir, "debug", "std::string", this);

	verbose_level = ParameterHolder::create<size_t>("verbose_level", "Verbose level [0-3]", verbose_level, 1, "size_t", this);
	testing_level = ParameterHolder::create<size_t>("testing_level", "Testing level [0-3]", testing_level, 0, "size_t", this);
	measure_level = ParameterHolder::create<size_t>("measure_level", "Measure level [0-3]", measure_level, 0, "size_t", this);

	print_matrices = ParameterHolder::create<bool>("print_matrices", "Print assembler matrices", print_matrices, false, "bool", this);
	remove_old_results = ParameterHolder::create<bool>("remove_old_results", "Keep only the last results", remove_old_results, false, "bool", this);

	environment = this;
}

}



