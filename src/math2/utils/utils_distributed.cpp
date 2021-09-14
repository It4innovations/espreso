
#include "utils_distributed.h"

#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "wrappers/mpi/communication.h"

#include <cstring>

using namespace espreso;

void DataSynchronization::gatherFromUpper(double *vals)
{
	for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
		memcpy(sBuffer[n].data(), vals + nStart[n], sizeof(double) * sBuffer[n].size());
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("receive MatrixCSRDistribution data.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		for (size_t i = 0; i < rIndices[n].size(); ++i) {
			vals[rIndices[n][i]] += rBuffer[n][i];
		}
	}
}

void DataSynchronization::scatterToUpper(double *vals)
{

}
