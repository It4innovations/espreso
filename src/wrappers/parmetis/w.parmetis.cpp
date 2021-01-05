
#ifdef HAVE_PARMETIS
#include "parmetis.h"
#endif

#include "w.parmetis.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

bool ParMETIS::islinked()
{
#ifdef HAVE_PARMETIS
	return true;
#endif
	return false;
}

esint ParMETIS::call(
			ParMETIS::METHOD method,
			MPIGroup &group,
			esint *edistribution,
			esint *eframes, esint *eneighbors,
			esint dimensions, float *coordinates,
			esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
			esint *partition)
{
	esint edgecut = 0;

#ifndef HAVE_PARMETIS
	eslog::globalerror("ESPRESO run-time error: cannot call ParMETIS library (the library is not linked).\n");
#else
	verticesWeightCount = std::max((esint)1, verticesWeightCount);

	esint wgtflag = 0;
	esint numflag = 0;
	esint parts = MPITools::procs->size;
	std::vector<float> partFraction(verticesWeightCount * parts, 1.0 / parts);
	std::vector<float> unbalanceTolerance(verticesWeightCount, 1.02);
	esint options[4] = { 0, 0, 0, PARMETIS_PSR_UNCOUPLED };
	float itr = 1e6;

	if (verticesWeights != NULL) {
		wgtflag += 2;
	}
	if (edgeWeights != NULL) {
		wgtflag += 1;
	}

	switch (method) {

	case ParMETIS::METHOD::ParMETIS_V3_PartKway:
		if (coordinates != NULL) {
			if (METIS_OK != ParMETIS_V3_PartGeomKway(
					edistribution,
					eframes, eneighbors,
					verticesWeights, edgeWeights,
					&wgtflag, &numflag, &dimensions, coordinates, &verticesWeightCount,
					&parts, partFraction.data(), unbalanceTolerance.data(),
					options,
					&edgecut, partition,
					&group.communicator)) {

				eslog::error("PARMETIS_ERROR while partitiate mesh to MPI processes by KWay utilizing coordinates.\n");
			}
		} else {
			if (METIS_OK != ParMETIS_V3_PartKway(
					edistribution,
					eframes, eneighbors,
					verticesWeights, edgeWeights,
					&wgtflag, &numflag, &verticesWeightCount,
					&parts, partFraction.data(), unbalanceTolerance.data(),
					options,
					&edgecut, partition,
					&group.communicator)) {

				eslog::error("PARMETIS_ERROR while partitiate mesh to MPI processes by KWay.\n");
			}
		}
		break;

	case ParMETIS::METHOD::ParMETIS_V3_RefineKway:
		if (METIS_OK != ParMETIS_V3_RefineKway(
				edistribution,
				eframes, eneighbors,
				verticesWeights, edgeWeights,
				&wgtflag, &numflag, &verticesWeightCount,
				&parts, partFraction.data(), unbalanceTolerance.data(),
				options,
				&edgecut, partition,
				&group.communicator)) {

			eslog::error("PARMETIS_ERROR while refine mesh partition to MPI processes by KWay.\n");
		}
		break;

	case ParMETIS::METHOD::ParMETIS_V3_AdaptiveRepart:
		if (METIS_OK != ParMETIS_V3_AdaptiveRepart(
				edistribution,
				eframes, eneighbors,
				verticesWeights, NULL, edgeWeights,
				&wgtflag, &numflag, &verticesWeightCount,
				&parts, partFraction.data(), unbalanceTolerance.data(), &itr,
				options,
				&edgecut, partition,
				&group.communicator)) {

			eslog::error("PARMETIS_ERROR while adaptive repartition mesh to MPI processes by KWay.\n");
		}
		break;

	case ParMETIS::METHOD::ParMETIS_V3_PartGeom:
		if (coordinates == NULL) {
			eslog::error("PARMETIS_ERROR:: PartGeom needs coordinates.\n");
		}
		if (METIS_OK != ParMETIS_V3_PartGeom(
				edistribution,
				&dimensions, coordinates,
				partition,
				&group.communicator)) {

			eslog::error("PARMETIS_ERROR while refine mesh partition to MPI processes by KWay.\n");
		}

	}
#endif

	return edgecut;
}


