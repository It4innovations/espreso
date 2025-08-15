
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPARSESOLVER_CPU_CONFIG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPARSESOLVER_CPU_CONFIG_H_

#include "config/description.h"

namespace espreso {

struct DualopTotalfetiImplicitGeneralSparseSolverCpuConfig: public ECFDescription {

    enum struct AUTOBOOL {
        AUTO,
        TRUE,
        FALSE
    };

    enum struct SPARSE_SOLVER_IMPL {
        AUTO,
        MKLPARDISO,
        SUITESPARSE,
        MUMPS,
        STRUMPACK,
        PASTIX,
        SUPERLU_DIST
    };

    AUTOBOOL parallel_set;
    AUTOBOOL parallel_update;
    AUTOBOOL parallel_apply;
    AUTOBOOL timers_outer;
    AUTOBOOL timers_inner;
    AUTOBOOL print_config;
    SPARSE_SOLVER_IMPL sparse_solver_impl;

    DualopTotalfetiImplicitGeneralSparseSolverCpuConfig();

};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPARSESOLVER_CPU_CONFIG_H_ */
