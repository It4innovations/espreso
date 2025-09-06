
#ifndef SRC_CONFIG_ECF_OPERATIONS_GPU_SCHUR_HCSX_DDNY_MANUAL_SIMPLE_H_
#define SRC_CONFIG_ECF_OPERATIONS_GPU_SCHUR_HCSX_DDNY_MANUAL_SIMPLE_H_

#include "config/description.h"

namespace espreso {

struct GpuSchurHcsxDdnyManualSimpleConfig: public ECFDescription {

    enum struct MATRIX_ORDER {
        AUTO,
        ROW_MAJOR,
        COL_MAJOR
    };

    enum struct SPDN {
        AUTO,
        SPARSE,
        DENSE
    };

    // make all available, throw the error that they dont provide factors later
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
    SPDN trsm_factor_spdn;
    SPARSE_SOLVER_IMPL sparse_solver_impl;

    GpuSchurHcsxDdnyManualSimpleConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_GPU_SCHUR_HCSX_DDNY_MANUAL_SIMPLE_H_ */
