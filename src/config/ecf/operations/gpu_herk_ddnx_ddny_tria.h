
#ifndef SRC_CONFIG_ECF_OPERATIONS_GPU_HERK_DDNX_DDNY_TRIA_H_
#define SRC_CONFIG_ECF_OPERATIONS_GPU_HERK_DDNX_DDNY_TRIA_H_

#include "config/description.h"

namespace espreso {

struct GpuHerkDdnxDdnyTriaConfig: public ECFDescription {

    enum struct HERK_TRIA_STRATEGY {
        AUTO,
        STAIRS,
        SQUARES
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

    struct GpuHerkDdnxDdnyTriaPartitionConfig: public ECFDescription {

        PARTITION_ALGORITHM algorithm;
        PARTITION_STRATEGY strategy;
        int chunk_size; // 0 = auto
        int chunk_count; // 0 = auto

        GpuHerkDdnxDdnyTriaPartitionConfig();
    };

    HERK_TRIA_STRATEGY strategy;
    GpuHerkDdnxDdnyTriaPartitionConfig partition;

    GpuHerkDdnxDdnyTriaConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_GPU_HERK_DDNX_DDNY_TRIA_H_ */
