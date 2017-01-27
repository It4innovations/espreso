
#ifndef SRC_CONFIG_ENVIRONMENT_H_
#define SRC_CONFIG_ENVIRONMENT_H_

#include "configuration.h"

namespace espreso {

struct Environment: public Configuration {

	Environment();

	int MPIrank = 0;
	int MPIsize = 1;

	std::string executable;

	PARAMETER(size_t, MKL_NUM_THREADS, "Number of MKL threads."                      , 1);
	PARAMETER(size_t, OMP_NUM_THREADS, "Number of OMP threads."                      , 1);
	PARAMETER(size_t, SOLVER_NUM_THREADS, "Number of threads used in ESPRESO solver.", 1);
	PARAMETER(size_t, PAR_NUM_THREADS, "Number of parallel threads."                 , 1);
	PARAMETER(size_t, CILK_NWORKERS, "Number of cilk++ threads."                     , 1);

	PARAMETER(std::string, log_dir, "Log directory.", "log");

	PARAMETER(size_t, verbose_level, "Verbose level [0-3]", 1);
	PARAMETER(size_t, testing_level, "Testing level [0-3]", 0);
	PARAMETER(size_t, measure_level, "Measure level [0-3]", 0);

	PARAMETER(bool, print_matrices, "Print assembler matrices.", false);
};

extern Environment *environment;

}



#endif /* SRC_CONFIG_ENVIRONMENT_H_ */
