
#include <signal.h>
#include <csignal>

#include "mpi.h"

#include "../configuration/globalconfiguration.h"
#include "factory/factory.h"

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

	ESINFO(OVERVIEW) << "Run ESPRESO on " << environment->MPIsize << " process(es).";

	Factory factory(configuration);

	factory.solve();
	factory.check(configuration.results);
	factory.finalize();

	MPI_Barrier(environment->MPICommunicator);
	MPI_Finalize();

	return 0;
}


