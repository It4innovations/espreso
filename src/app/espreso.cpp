
#include <signal.h>
#include <csignal>

#include "mpi.h"

#include "../configuration/globalconfiguration.h"
#include "factory/factory.h"

#ifdef READEX_LEVEL_1
#include <readex.h>
#include <readex_regions.h>
#endif
using namespace espreso;

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGSEGV:
		ESINFO(ERROR) << "Invalid memory reference";
		break;
	case SIGFPE:
		ESINFO(ERROR) << "Erroneous arithmetic operation";
		break;
	}
}


int main(int argc, char **argv)
{
	std::signal(SIGFPE, signalHandler);
	std::signal(SIGSEGV, signalHandler);

	MPI_Init(&argc, &argv);

#ifdef READEX_LEVEL_1
	READEX_INIT();
	READEX_PHASE_START(REG_Main, "Main", SCOREP_USER_REGION_TYPE_PHASE);
#endif

	GlobalConfiguration configuration(&argc, &argv);

	ESINFO(OVERVIEW) << "Run ESPRESO on " << environment->MPIsize << " process(es).";

	Factory factory(configuration);

	factory.solve();
	factory.finalize();

	MPI_Barrier(MPI_COMM_WORLD);
#ifdef READEX_LEVEL_1
	READEX_PHASE_STOP(REG_Main);
	READEX_CLOSE();
#endif
	MPI_Finalize();

	return 0;
}


