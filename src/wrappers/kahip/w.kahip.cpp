
#ifdef HAVE_KAHIP
#include "kaHIP_interface.h"
#endif

#include "w.kahip.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

bool KaHIP::islinked()
{
#ifdef HAVE_KAHIP
	return true;
#endif
	return false;
}

esint KaHIP::call(
		const KaHIPConfiguration &options,
		esint verticesCount,
		esint *eframes, esint *eneighbors,
		esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
		esint parts, esint *partition)
{
	esint edgecut = 0;

#ifndef HAVE_KAHIP
	eslog::globalerror("ESPRESO run-time error: cannot call KaHIP library (the library is not linked).\n");
#else
	double imbalance = 0.03;

	kaffpa(&verticesCount, verticesWeights, eframes, edgeWeights, eneighbors, &parts, &imbalance, true, 0, STRONG, &edgecut, partition);
#endif
	return edgecut;
}


