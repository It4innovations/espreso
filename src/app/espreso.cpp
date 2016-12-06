
#include <signal.h>
#include <csignal>

#include "../config/description.h"
#include "factory/factory.h"
#include "esinput.h"


using namespace espreso;

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGSEGV:
		ESINFO(ERROR) << "Invalid memory reference";
		break;
	case SIGFPE:
		ESINFO(ERROR) << "Errorneous arithmetic operation";
		break;
	}
}


int main(int argc, char **argv)
{
	std::signal(SIGFPE, signalHandler);
	std::signal(SIGSEGV, signalHandler);

	MPI_Init(&argc, &argv);

	GlobalConfiguration configuration(&argc, &argv);
	ParametersReader::fromArguments(&argc, &argv);

	ESINFO(OVERVIEW) << "Run ESPRESO on " << config::env::MPIsize << " process(es).";
	ParametersReader::printParameters(config::parameters, config::info::VERBOSE_LEVEL);

	Factory factory(configuration);
	factory.solve("result");

	ESTEST(EVALUATION)
		<< (fabs(factory.norm() - config::solver::NORM) > 1e-3 && !config::env::MPIrank ? TEST_FAILED : TEST_PASSED)
		<< "Norm of the solution " << factory.norm() << " is not " << config::solver::NORM << ".";

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}


