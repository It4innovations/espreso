
#include "envinfo.h"

#include <sstream>
#include <cstdlib>
#include <omp.h>

int espreso::info::env::MKL_NUM_THREADS = 1;
int espreso::info::env::OMP_NUM_THREADS = 1;
int espreso::info::env::SOLVER_NUM_THREADS = 1;
int espreso::info::env::PAR_NUM_THREADS = 1;
int espreso::info::env::threads = 1;

void espreso::info::env::set()
{
	auto getEnv = [] (int &value, const char *name)
	{
		char *var = getenv(name);
		if (var != NULL) {
			std::stringstream ss(var);
			ss >> value;
		}
	};

	getEnv(MKL_NUM_THREADS   , "MKL_NUM_THREADS");
	getEnv(OMP_NUM_THREADS   , "OMP_NUM_THREADS");
	getEnv(SOLVER_NUM_THREADS, "SOLVER_NUM_THREADS");
	getEnv(PAR_NUM_THREADS   , "PAR_NUM_THREADS");

	omp_set_num_threads(OMP_NUM_THREADS);
	threads = OMP_NUM_THREADS;
}

char* espreso::info::env::pwd()
{
	return getenv("PWD");
}


