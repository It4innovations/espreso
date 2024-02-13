
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_EXPLICIT_GPU_CONFIG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_EXPLICIT_GPU_CONFIG_H_

#include "config/description.h"

namespace espreso {

struct DualOperatorExplicitGpuConfig: public ECFDescription {

    enum class CONCURRENCY {
        PARALLEL,
        SEQ_CONTINUE,
        SEQ_WAIT
    };

    enum class MATRIX_STORAGE {
        CSR,
        DENSE_ROW_MAJOR
    };

    enum class TRSM1_SOLVE_TYPE {
        L,
        LHH
    };

    enum class TRSM2_SOLVE_TYPE {
        U,
        UHH
    };

    enum class MATRIX_ORDER {
        ROW_MAJOR,
        COL_MAJOR
    };

    enum class PATH_IF_HERMITIAN {
        SECOND_TRSM,
        HERK
    };

    enum class QUEUE_COUNT {
        PER_THREAD,
        PER_DOMAIN
    };

    enum class DEVICE {
        CPU,
        GPU
    };

    CONCURRENCY concurrency_set;
    CONCURRENCY concurrency_update;
    CONCURRENCY concurrency_apply;
    MATRIX_STORAGE trsm1_factor_storage;
    MATRIX_STORAGE trsm2_factor_storage;
    TRSM1_SOLVE_TYPE trsm1_solve_type;
    TRSM2_SOLVE_TYPE trsm2_solve_type;
    MATRIX_ORDER trsm_rhs_sol_order;
    PATH_IF_HERMITIAN path_if_hermitian;
    QUEUE_COUNT queue_count;
    DEVICE apply_scatter_gather_where;

	DualOperatorExplicitGpuConfig();

};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_EXPLICIT_GPU_CONFIG_H_ */
