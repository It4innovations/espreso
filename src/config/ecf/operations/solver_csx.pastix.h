
#ifndef SRC_CONFIG_ECF_OPERATIONS_SOLVER_CSX_PASTIX_H_
#define SRC_CONFIG_ECF_OPERATIONS_SOLVER_CSX_PASTIX_H_

#include "config/description.h"

namespace espreso {

struct SolverCsxPastixConfig: public ECFDescription {

    enum struct AUTOBOOL {
        AUTO,
        TRUE,
        FALSE
    };

    AUTOBOOL use_gpu;

    SolverCsxPastixConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_SOLVER_CSX_PASTIX_H_ */
