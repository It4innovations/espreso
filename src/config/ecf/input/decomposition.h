
#ifndef SRC_CONFIG_ECF_INPUT_DECOMPOSITION_H_
#define SRC_CONFIG_ECF_INPUT_DECOMPOSITION_H_

#include "config/description.h"

#include <string>

namespace espreso {

struct ParMETISConfiguration: public ECFDescription {

    bool refinement;
    double tolerance;

    ParMETISConfiguration();
};

struct METISConfiguration: public ECFDescription {

    enum class PTYPE {
        RB, KWAY
    };

    enum class OBJTYPE {
        CUT, VOL, NODE
    };

    enum class CTYPE {
        RM, SHEM
    };

    enum class IPTYPE {
        GROW, RANDOM, EDGE, NODE
    };

    enum class RTYPE {
        FM, GREEDY, SEP2SIDED, SEP1SIDED
    };

    PTYPE ptype;
    OBJTYPE objtype;
    CTYPE ctype;
    IPTYPE iptype;
    RTYPE rtype;

    int ncuts, nseps, niter, seed, minconn, no2hp, contig, compress, ccorder, dbglvl;
    double pfactor, ufactor;

    METISConfiguration();
};

struct PTScotchConfiguration: public ECFDescription {

    PTScotchConfiguration();
};

struct ScotchConfiguration: public ECFDescription {

    ScotchConfiguration();
};

struct KaHIPConfiguration: public ECFDescription {

    KaHIPConfiguration();
};

struct DecompositionConfiguration: public ECFDescription {

    enum class ParallelDecomposer {
        NONE,
        METIS,
        PARMETIS,
        PTSCOTCH,
        HILBERT_CURVE
    };

    enum class SequentialDecomposer {
        NONE,
        METIS,
        SCOTCH,
        KAHIP
    };

    ParallelDecomposer parallel_decomposer;
    SequentialDecomposer sequential_decomposer;
    int mesh_duplication;
    int domains;

    bool force_continuity;
    bool separate_materials, separate_regions, separate_etypes;
    ParMETISConfiguration parmetis_options;
    METISConfiguration metis_options;
    PTScotchConfiguration ptscotch_options;
    ScotchConfiguration scotch_options;
    KaHIPConfiguration kahip_options;

    DecompositionConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_DECOMPOSITION_H_ */
