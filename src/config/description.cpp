
#include "description.h"

namespace espreso {

Environment *environment;
OutputConfiguration *output;

Environment::Environment(): executable("espreso")
{
	int initialized;
	MPI_Initialized(&initialized);
	if (initialized) {
		MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
		MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	}
	environment = this;
}

OutputConfiguration::OutputConfiguration()
{
	output = this;
}

}



