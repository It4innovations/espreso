
#ifdef HAVE_METIS
#include "metis.h"
#endif

#include "w.metis.h"
#include "config/ecf/input/decomposition.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

bool METIS::islinked()
{
#ifdef HAVE_METIS
	return true;
#else
	return false;
#endif
}

esint METIS::call(
		const METISConfiguration &options,
		esint verticesCount,
		esint *eframes, esint *eneighbors,
		esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
		esint parts, esint *partition)
{
	esint edgecut = 0;

#ifndef HAVE_METIS
	eslog::globalerror("ESPRESO run-time error: cannot call METIS library (the library is not linked).\n");
#else
	verticesWeightCount = std::max((esint)1, verticesWeightCount);

	esint moptions[METIS_NOPTIONS];
	METIS_SetDefaultOptions(moptions);

	// HEURISTICS
	switch (options.objective_type) {
	case METISConfiguration::OBJECTIVE_TYPE::VOLUME:
		moptions[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_VOL;
		break;
	case METISConfiguration::OBJECTIVE_TYPE::EDGECUT:
		moptions[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
		break;
	default:
		eslog::error("METIS WRAPPER error: invalid objective type.\n");
	}

	moptions[METIS_OPTION_CTYPE]     = METIS_CTYPE_SHEM;
	moptions[METIS_OPTION_IPTYPE]    = METIS_IPTYPE_GROW;
	moptions[METIS_OPTION_RTYPE]     = METIS_RTYPE_FM;
	moptions[METIS_OPTION_NITER]     = 20;
	moptions[METIS_OPTION_UFACTOR]   = 100; // imbalance (1 + x) / 1000

	// MANDATORY
	moptions[METIS_OPTION_CONTIG]    = options.continuous;
	moptions[METIS_OPTION_MINCONN]   = 0;
	moptions[METIS_OPTION_NUMBERING] = 0;
	moptions[METIS_OPTION_DBGLVL]    = 0;

	// KEPT DEFAULT
	// options[METIS_OPTION_NO2HOP]    = 0;
	// options[METIS_OPTION_NCUTS]     = 1;
	// options[METIS_OPTION_SEED]      = 0;

	if (parts > 1) {
		if (METIS_OK != METIS_PartGraphKway(
				&verticesCount, &verticesWeightCount,
				eframes, eneighbors,
				verticesWeights, NULL, edgeWeights,
				&parts, NULL, NULL,
				moptions, &edgecut, partition)) {

			eslog::error("METIS_ERROR while KWay decomposition.\n");
		}
	} else {
		std::fill(partition, partition + verticesCount, 0);
	}

#endif

	return edgecut;
}


