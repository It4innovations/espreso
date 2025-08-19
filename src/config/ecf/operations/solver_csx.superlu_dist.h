
#ifndef SRC_CONFIG_ECF_OPERATIONS_SOLVER_CSX_SUPERLU_DIST_H_
#define SRC_CONFIG_ECF_OPERATIONS_SOLVER_CSX_SUPERLU_DIST_H_

#include "config/description.h"

namespace espreso {

struct SolverCsxSuperludistConfig: public ECFDescription {

    enum struct AUTOBOOL {
        AUTO,
        TRUE,
        FALSE
    };

    AUTOBOOL use_gpu;

    SolverCsxSuperludistConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_SOLVER_CSX_SUPERLU_DIST_H_ */
