
#include "system.h"
#include "config/configuration.hpp"

using namespace espreso;

SystemConfiguration::SystemConfiguration()
{
	mpi_affinity = 0;

	REGISTER(mpi_affinity, ECFMetaData()
			.setdescription({ "Number of hwthreads assigned to each MPI process." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));
}


