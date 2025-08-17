
#include "config/ecf/operations/gpu_trsm_hcsx_ddny_tria.h"

#include "config/configuration.hpp"

using namespace espreso;

GpuTrsmHcsxDdnyTriaConfig::GpuTrsmHcsxDdnyTriaConfig()
{
    strategy = TRSM_TRIA_STRATEGY::AUTO;
    REGISTER(strategy, ECFMetaData()
        .setdescription({ "Splitting strategy" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("SPLIT_RHS").setdescription("Split RHS/SOL matrix into block columns"))
        .addoption(ECFOption().setname("SPLIT_FACTOR").setdescription("Split factor into blocks"))
    );

    REGISTER(partition, ECFMetaData()
            .setdescription({ "partition config" }));

    REGISTER(split_rhs_config, ECFMetaData()
            .setdescription({ "config for the split_rhs strategy" }));

    REGISTER(split_factor_config, ECFMetaData()
            .setdescription({ "config for the split_factor strategy" }));
}

GpuTrsmHcsxDdnyTriaConfig::GpuTrsmHcsxDdnyTriaPartitionConfig::GpuTrsmHcsxDdnyTriaPartitionConfig()
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

GpuTrsmHcsxDdnyTriaConfig::GpuTrsmHcsxDdnyTriaSplitrhsConfig::GpuTrsmHcsxDdnyTriaSplitrhsConfig()
{
    factor_order_sp = MATRIX_ORDER::AUTO;
    REGISTER(factor_order_sp, ECFMetaData()
        .setdescription({ "Memory order of sparse factor part" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    factor_order_dn = MATRIX_ORDER::AUTO;
    REGISTER(factor_order_dn, ECFMetaData()
        .setdescription({ "Memory order of dense factor part" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    spdn_criteria = SPDN_CRITERIA::AUTO;
    REGISTER(spdn_criteria, ECFMetaData()
        .setdescription({ "Criterea for choosing between sparse and dense chunk of factor" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("SPARSE_ONLY").setdescription("Use only sparse factor"))
        .addoption(ECFOption().setname("DENSE_ONLY").setdescription("Use only dense factor"))
        .addoption(ECFOption().setname("FRACTION_OF_NUM_CHUNKS").setdescription("Given fraction of chunks will be sparse"))
        .addoption(ECFOption().setname("FRACTION_OF_FACTOR_SIZE").setdescription("Given fraction of factor size will be sparse"))
        .addoption(ECFOption().setname("FACTOR_DENSITY").setdescription("Chunks with higher density will be dense"))
    );

    spdn_param_frac_of_num_chunks = 0;
    REGISTER(spdn_param_frac_of_num_chunks, ECFMetaData()
            .setdescription({ "Parameter for the FRACTION_OF_NUM_CHUNKS sparse-dense criteria. (0.0-1.0)" })
            .setdatatype({ ECFDataType::FLOAT }));

    spdn_param_frac_of_factor_size = 0;
    REGISTER(spdn_param_frac_of_factor_size, ECFMetaData()
            .setdescription({ "Parameter for the FRACTION_OF_FACTOR_SIZE sparse-dense criteria. (0.0-1.0)" })
            .setdatatype({ ECFDataType::FLOAT }));

    spdn_param_factor_density = 0;
    REGISTER(spdn_param_factor_density, ECFMetaData()
            .setdescription({ "Parameter for the FACTOR_DENSITY sparse-dense criteria. (0.0-1.0)" })
            .setdatatype({ ECFDataType::FLOAT }));
}

GpuTrsmHcsxDdnyTriaConfig::GpuTrsmHcsxDdnyTriaSplitfactorConfig::GpuTrsmHcsxDdnyTriaSplitfactorConfig()
{
    trsm_factor_spdn = SPDN::AUTO;
    REGISTER(trsm_factor_spdn, ECFMetaData()
        .setdescription({ "TRSM kernel will use sparse or dense factor submatrix" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("SPARSE").setdescription("Sparse factor submatrix"))
        .addoption(ECFOption().setname("DENSE").setdescription("Dense factor submatrix"))
    );

    trsm_factor_order = MATRIX_ORDER::AUTO;
    REGISTER(trsm_factor_order, ECFMetaData()
        .setdescription({ "Memory order of factor submatrix for the TRSM kernel" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    gemm_factor_pruning = PRUNING_STRATEGY::AUTO;
    REGISTER(gemm_factor_pruning, ECFMetaData()
        .setdescription({ "Pruning used in the GEMM kernel (disregarding empty rows/cols)" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("NO_PRUNING").setdescription("Do not prune at all"))
        .addoption(ECFOption().setname("ROWS_ONLY").setdescription("Prune only empty rows"))
        .addoption(ECFOption().setname("COLS_ONLY").setdescription("Prune only empty columns"))
        .addoption(ECFOption().setname("ROWS_AND_COLS").setdescription("Prune both rows and columns"))
    );

    gemm_factor_order_sp = MATRIX_ORDER::AUTO;
    REGISTER(gemm_factor_order_sp, ECFMetaData()
        .setdescription({ "Memory order of sparse factor submatrix for the GEMM kernel" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    gemm_factor_order_dn = MATRIX_ORDER::AUTO;
    REGISTER(gemm_factor_order_dn, ECFMetaData()
        .setdescription({ "Memory order of dense factor submatrix for the GEMM kernel" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    gemm_spdn_criteria = SPDN_CRITERIA::AUTO;
    REGISTER(gemm_spdn_criteria, ECFMetaData()
        .setdescription({ "Criteria for decision between sparse and dense factor submatrix in GEMM kernel" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("SPARSE_ONLY").setdescription("Use only sparse factor"))
        .addoption(ECFOption().setname("DENSE_ONLY").setdescription("Use only dense factor"))
        .addoption(ECFOption().setname("FRACTION_OF_NUM_CHUNKS").setdescription("Given fraction of chunks (from the left) will use sparse"))
        .addoption(ECFOption().setname("FRACTION_OF_FACTOR_SIZE").setdescription("Given fraction of factor size (from the left) will use sparse"))
        .addoption(ECFOption().setname("FACTOR_DENSITY").setdescription("Chunks with higher density will be dense"))
    );

    spdn_param_frac_of_num_chunks = 0;
    REGISTER(spdn_param_frac_of_num_chunks, ECFMetaData()
            .setdescription({ "Parameter for the FRACTION_OF_NUM_CHUNKS sparse-dense criteria. (0.0-1.0)" })
            .setdatatype({ ECFDataType::FLOAT }));

    spdn_param_frac_of_factor_size = 0;
    REGISTER(spdn_param_frac_of_factor_size, ECFMetaData()
            .setdescription({ "Parameter for the FRACTION_OF_FACTOR_SIZE sparse-dense criteria. (0.0-1.0)" })
            .setdatatype({ ECFDataType::FLOAT }));

    spdn_param_factor_density = 0;
    REGISTER(spdn_param_factor_density, ECFMetaData()
            .setdescription({ "Parameter for the FACTOR_DENSITY sparse-dense criteria. (0.0-1.0)" })
            .setdatatype({ ECFDataType::FLOAT }));
}
