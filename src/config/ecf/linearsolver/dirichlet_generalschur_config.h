
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_DIRICHLET_GENERALSCHUR_CONFIG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_DIRICHLET_GENERALSCHUR_CONFIG_H_

#include "config/description.h"

namespace espreso {

struct DirichletGeneralSchurConfig: public ECFDescription {

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

    enum struct SCHUR_IMPL_CPU {
        AUTO,
        MANUAL_SIMPLE,
        TRIANGULAR,
        MKLPARDISO,
        SPARSE_SOLVER,
        MUMPS,
        PASTIX,
    };

    enum struct SCHUR_IMPL_GPU {
        AUTO,
        MANUAL_SIMPLE,
        TRIANGULAR,
    };

    AUTOBOOL parallel_set;
    AUTOBOOL parallel_update;
    AUTOBOOL parallel_apply;
    AUTOBOOL timers_outer;
    AUTOBOOL timers_inner;
    AUTOBOOL print_config;
    MATRIX_ORDER order_sc;
    SCHUR_IMPL_CPU schur_impl_cpu;
    SCHUR_IMPL_GPU schur_impl_gpu;

    DirichletGeneralSchurConfig();

};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_DIRICHLET_GENERALSCHUR_CONFIG_H_ */
