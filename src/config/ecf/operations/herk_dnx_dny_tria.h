
#ifndef SRC_CONFIG_ECF_OPERATIONS_HERK_DNX_DNY_TRIA_H_
#define SRC_CONFIG_ECF_OPERATIONS_HERK_DNX_DNY_TRIA_H_

#include "config/description.h"

namespace espreso {

struct HerkDnxDnyTriaConfig: public ECFDescription {

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

    struct HerkDnxDnyTriaPartitionConfig: public ECFDescription {

        PARTITION_ALGORITHM algorithm;
        PARTITION_STRATEGY strategy;
        int chunk_size; // 0 = auto
        int chunk_count; // 0 = auto

        HerkDnxDnyTriaPartitionConfig();
    };

    HERK_TRIA_STRATEGY strategy;
    HerkDnxDnyTriaPartitionConfig partition;

    HerkDnxDnyTriaConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_HERK_DNX_DNY_TRIA_H_ */
