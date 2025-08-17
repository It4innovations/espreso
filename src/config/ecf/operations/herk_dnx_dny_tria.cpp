
#include "config/ecf/operations/herk_dnx_dny_tria.h"

#include "config/configuration.hpp"

using namespace espreso;

HerkDnxDnyTriaConfig::HerkDnxDnyTriaConfig()
{
    strategy = HERK_TRIA_STRATEGY::AUTO;
    REGISTER(strategy, ECFMetaData()
        .setdescription({ "Splitting strategy" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("STAIRS").setdescription("Split input matrix in the n direction, output resembles stairs"))
        .addoption(ECFOption().setname("SQUARES").setdescription("Split input matrix in the k direction, output resembles squares"))
    );

    REGISTER(partition, ECFMetaData()
            .setdescription({ "partition config" }));
}

HerkDnxDnyTriaConfig::HerkDnxDnyTriaPartitionConfig::HerkDnxDnyTriaPartitionConfig()
{
    algorithm = PARTITION_ALGORITHM::AUTO;
    REGISTER(algorithm, ECFMetaData()
        .setdescription({ "Partition algorithm" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("UNIFORM").setdescription("Chunks will have uniform size"))
        .addoption(ECFOption().setname("MINIMUM_WORK").setdescription("Chunk bounds are determined to minimize FLOP count"))
    );

    strategy = PARTITION_STRATEGY::AUTO;
    REGISTER(strategy, ECFMetaData()
        .setdescription({ "Partition strategy (chunk scaling)" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("CHUNK_SIZE").setdescription("Chunk size is constant"))
        .addoption(ECFOption().setname("CHUNK_COUNT").setdescription("Chunk count is constant"))
    );

    chunk_size = 0;
    REGISTER(chunk_size, ECFMetaData()
            .setdescription({ "Size of chunk for the CHUNK_SIZE partition strategy" })
            .setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

    chunk_count = 0;
    REGISTER(chunk_count, ECFMetaData()
            .setdescription({ "Number of chunks for the CHUNK_COUNT partition strategy" })
            .setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));
}
