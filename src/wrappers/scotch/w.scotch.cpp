
#include "w.scotch.h"
#include "basis/utilities/communication.h"
#include "esinfo/eslog.hpp"

#ifdef HAVE_SCOTCH
#include "ptscotch.h"
#endif

using namespace espreso;

bool Scotch::islinked()
{
#ifdef HAVE_SCOTCH
	return true;
#endif
	return false;
}

esint Scotch::call(
		const ScotchConfiguration &options,
		esint verticesCount,
		esint *eframes, esint *eneighbors,
		esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
		esint parts, esint *partition)
{
	esint edgecut = 0;

#ifndef HAVE_SCOTCH
	eslog::globalerror("ESPRESO run-time error: cannot call PT-Scotch library (the library is not linked).\n");
#else
	SCOTCH_Dgraph grafdat;
	if (SCOTCH_dgraphInit(&grafdat, MPI_COMM_SELF)) {
		eslog::error("PTSCOTCH_ERROR while call 'SCOTCH_dgraphInit'.\n");
	}

	esint baseval = 0;
	esint vertlocnbr = verticesCount;
	esint vertlocmax = vertlocnbr;
	esint *vertloctab = eframes;
	esint *vendloctab = NULL;
	esint *veloloctab = NULL;
	esint *vlblloctab = NULL;
	esint edgelocnbr = eframes[vertlocnbr];
	esint edgelocsiz = edgelocnbr;
	esint *edgeloctab = eneighbors;
	esint *edgegsttab = NULL;
	esint *edloloctab = NULL;
	int ierr;

	if ((ierr = SCOTCH_dgraphBuild(
			&grafdat, baseval,
			vertlocnbr, vertlocmax, vertloctab, vendloctab, veloloctab, vlblloctab,
			edgelocnbr, edgelocsiz, edgeloctab, edgegsttab, edloloctab))) {

		eslog::error("PTSCOTCH_ERROR while call 'SCOTCH_dgraphBuild() = %d'.\n", ierr);
	}

	SCOTCH_Strat straptr;

	if (SCOTCH_stratInit(&straptr)) {
		eslog::error("PTSCOTCH_ERROR while call 'SCOTCH_stratInit() = %d'.\n", ierr);
	}

	if ((ierr = SCOTCH_dgraphPart(&grafdat, parts, &straptr, partition))) {
		eslog::error("PTSCOTCH_ERROR while call 'SCOTCH_dgraphPart() = %d'.\n", ierr);
	}

	SCOTCH_dgraphExit(&grafdat);

#endif
	return edgecut;
}


