
#ifndef SRC_CONFIG_ECF_OPERATIONS_SCHUR_CSX_DNY_MANUAL_SIMPLE_H_
#define SRC_CONFIG_ECF_OPERATIONS_SCHUR_CSX_DNY_MANUAL_SIMPLE_H_

#include "config/description.h"

namespace espreso {

struct SchurCsxDnyManualSimpleConfig: public ECFDescription {

    enum struct MATRIX_ORDER {
        AUTO,
        ROW_MAJOR,
        COL_MAJOR
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

    MATRIX_ORDER order_X;
    MATRIX_ORDER order_L;
    SPARSE_SOLVER_IMPL sparse_solver_impl;

    SchurCsxDnyManualSimpleConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_SCHUR_CSX_DNY_MANUAL_SIMPLE_H_ */
