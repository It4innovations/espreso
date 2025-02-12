
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRE_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRE_H_

#include "hypreboomeramg.h"
#include "hyprepcg.h"
#include "hypregmres.h"
#include "hypreflexgmres.h"
#include "hyprelgmres.h"
#include "hyprebicgstab.h"
#include "hyprecgnr.h"

namespace espreso {

struct HYPREConfiguration: public ECFDescription {

    enum class SOLVER_TYPE {
        BoomerAMG,
        PCG,
        GMRES,
        FlexGMRES,
        LGMRES,
        BiCGSTAB,
        CGNR
    };

    SOLVER_TYPE solver_type;

    HYPREBoomerAMGConfiguration boomeramg;
    HYPREPCGConfiguration pcg;
    HYPREGMRESConfiguration gmres;
    HYPREFlexGMRESConfiguration flexgmres;
    HYPRELGMRESConfiguration lgmres;
    HYPREBiCGSTABConfiguration bicgstab;
    HYPRECGNRConfiguration cgnr;

    HYPREConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRE_H_ */
