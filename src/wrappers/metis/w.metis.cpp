
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

    switch (options.ptype) {
    case METISConfiguration::PTYPE::RB:   moptions[METIS_OPTION_PTYPE] = METIS_PTYPE_RB; break;
    case METISConfiguration::PTYPE::KWAY: moptions[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; break;
    }

    switch (options.objtype) {
    case METISConfiguration::OBJTYPE::CUT:  moptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; break;
    case METISConfiguration::OBJTYPE::VOL:  moptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; break;
    case METISConfiguration::OBJTYPE::NODE: moptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_NODE; break;
    }

    switch (options.ctype) {
    case METISConfiguration::CTYPE::RM:   moptions[METIS_OPTION_CTYPE] = METIS_CTYPE_RM; break;
    case METISConfiguration::CTYPE::SHEM: moptions[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; break;
    }

    switch (options.iptype) {
    case METISConfiguration::IPTYPE::GROW:   moptions[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW; break;
    case METISConfiguration::IPTYPE::RANDOM: moptions[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; break;
    case METISConfiguration::IPTYPE::EDGE:   moptions[METIS_OPTION_IPTYPE] = METIS_IPTYPE_EDGE; break;
    case METISConfiguration::IPTYPE::NODE:   moptions[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE; break;
    }

    switch (options.rtype) {
    case METISConfiguration::RTYPE::FM:        moptions[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; break;
    case METISConfiguration::RTYPE::GREEDY:    moptions[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY; break;
    case METISConfiguration::RTYPE::SEP2SIDED: moptions[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED; break;
    case METISConfiguration::RTYPE::SEP1SIDED: moptions[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP1SIDED; break;
    }

    moptions[METIS_OPTION_NCUTS]     = options.ncuts;
    moptions[METIS_OPTION_NITER]     = options.niter;
    moptions[METIS_OPTION_SEED]      = options.seed;
    moptions[METIS_OPTION_MINCONN]   = options.minconn;
    moptions[METIS_OPTION_NO2HOP]    = options.no2hp;
    moptions[METIS_OPTION_CONTIG]    = options.contig;
    moptions[METIS_OPTION_COMPRESS]  = options.compress;
    moptions[METIS_OPTION_CCORDER]   = options.ccorder;
    moptions[METIS_OPTION_PFACTOR]   = options.pfactor;
    moptions[METIS_OPTION_UFACTOR]   = options.ufactor;
    moptions[METIS_OPTION_DBGLVL]    = options.dbglvl;

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


