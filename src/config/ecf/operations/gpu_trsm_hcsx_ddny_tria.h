
#ifndef SRC_CONFIG_ECF_OPERATIONS_GPU_TRSM_HCSX_DDNY_TRIA_H_
#define SRC_CONFIG_ECF_OPERATIONS_GPU_TRSM_HCSX_DDNY_TRIA_H_

#include "config/description.h"

namespace espreso {

struct GpuTrsmHcsxDdnyTriaConfig: public ECFDescription {

    enum struct TRSM_TRIA_STRATEGY {
        AUTO,
        SPLIT_RHS,
        SPLIT_FACTOR
    };

    enum struct MATRIX_ORDER {
        AUTO,
        ROW_MAJOR,
        COL_MAJOR
    };

    enum struct SPDN_CRITERIA {
        AUTO,
        SPARSE_ONLY,
        DENSE_ONLY,
        FRACTION_OF_NUM_CHUNKS,
        FRACTION_OF_FACTOR_SIZE,
        FACTOR_DENSITY
    };

    enum struct PARTITION_ALGORITHM {
        AUTO,
        UNIFORM,
        MINIMUM_WORK
    };

    enum struct PARTITION_STRATEGY {
        AUTO,
        CHUNK_SIZE,
        CHUNK_COUNT
    };

    enum struct SPDN {
        AUTO,
        SPARSE,
        DENSE
    };

    enum struct PRUNING_STRATEGY {
        AUTO,
        NO_PRUNING,
        ROWS_ONLY,
        COLS_ONLY,
        ROWS_AND_COLS
    };

    struct GpuTrsmHcsxDdnyTriaPartitionConfig: public ECFDescription {

        PARTITION_ALGORITHM algorithm;
        PARTITION_STRATEGY strategy;
        int chunk_size; // 0 = auto
        int chunk_count; // 0 = auto

        GpuTrsmHcsxDdnyTriaPartitionConfig();
    };

    struct GpuTrsmHcsxDdnyTriaSplitrhsConfig: public ECFDescription {
        MATRIX_ORDER factor_order_sp;
        MATRIX_ORDER factor_order_dn;
        SPDN_CRITERIA spdn_criteria;
        double spdn_param_frac_of_num_chunks; // 0 = auto
        double spdn_param_frac_of_factor_size; // 0 = auto
        double spdn_param_factor_density; // 0 = auto

        GpuTrsmHcsxDdnyTriaSplitrhsConfig();
    };

    struct GpuTrsmHcsxDdnyTriaSplitfactorConfig: public ECFDescription {
        SPDN trsm_factor_spdn;
        MATRIX_ORDER trsm_factor_order;
        PRUNING_STRATEGY gemm_factor_pruning;
        MATRIX_ORDER gemm_factor_order_sp;
        MATRIX_ORDER gemm_factor_order_dn;
        SPDN_CRITERIA gemm_spdn_criteria;
        double spdn_param_frac_of_num_chunks; // 0 = auto
        double spdn_param_frac_of_factor_size; // 0 = auto
        double spdn_param_factor_density; // 0 = auto

        GpuTrsmHcsxDdnyTriaSplitfactorConfig();
    };

    TRSM_TRIA_STRATEGY strategy;
    GpuTrsmHcsxDdnyTriaPartitionConfig partition;
    GpuTrsmHcsxDdnyTriaSplitrhsConfig split_rhs_config;
    GpuTrsmHcsxDdnyTriaSplitfactorConfig split_factor_config;

    GpuTrsmHcsxDdnyTriaConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_GPU_TRSM_HCSX_DDNY_TRIA_H_ */
