
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERALSCHUR_GPU_CONFIG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERALSCHUR_GPU_CONFIG_H_

#include "config/description.h"

namespace espreso {

struct DualopTotalfetiExplicitGeneralSchurGpuConfig: public ECFDescription {

    enum struct AUTOBOOL {
        AUTO,
        TRUE,
        FALSE
    };

    enum struct MATRIX_ORDER {
        AUTO,
        ROW_MAJOR,
        COL_MAJOR
    };

    enum struct MAINLOOP_UPDATE_SPLIT {
        AUTO,
        COMBINED,
        SEPARATE
    };

    enum struct SCHUR_IMPL {
        AUTO,
        MANUAL_SIMPLE,
        TRIANGULAR
    };

    enum struct CPU_GPU {
        AUTO,
        CPU,
        GPU
    };

    AUTOBOOL parallel_set;
    AUTOBOOL parallel_update;
    AUTOBOOL parallel_apply;
    MAINLOOP_UPDATE_SPLIT mainloop_update_split;
    AUTOBOOL synchronize_after_update_mainloop;
    AUTOBOOL timers_outer;
    AUTOBOOL timers_inner;
    AUTOBOOL print_config;
    MATRIX_ORDER order_F;
    SCHUR_IMPL schur_impl;
    CPU_GPU apply_where;

    DualopTotalfetiExplicitGeneralSchurGpuConfig();

};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERALSCHUR_GPU_CONFIG_H_ */
