
#include "config/ecf/operations/gpu_schur_hcsx_ddny.manual_simple.h"

#include "config/configuration.hpp"

using namespace espreso;

GpuSchurHcsxDdnyManualSimpleConfig::GpuSchurHcsxDdnyManualSimpleConfig()
{
    order_X = MATRIX_ORDER::AUTO;
    REGISTER(order_X, ECFMetaData()
        .setdescription({ "Memory order of the intermediate dense matrix X" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    order_L = MATRIX_ORDER::AUTO;
    REGISTER(order_L, ECFMetaData()
        .setdescription({ "Memory order of the sparse factor L" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    trsm_factor_spdn = SPDN::AUTO;
    REGISTER(trsm_factor_spdn, ECFMetaData()
        .setdescription({ "Sparse or dense factor is used in TRSM kernel" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("SPARSE").setdescription("Sparse factor"))
        .addoption(ECFOption().setname("DENSE").setdescription("Dense factor"))
    );

    sparse_solver_impl = SPARSE_SOLVER_IMPL::AUTO;
    REGISTER(sparse_solver_impl, ECFMetaData()
        .setdescription({ "Sparse solver wrapper to be used" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("MKLPARDISO").setdescription("MKL Pardiso"))
        .addoption(ECFOption().setname("SUITESPARSE").setdescription("SuiteSparse"))
        .addoption(ECFOption().setname("MUMPS").setdescription("MUMPS"))
        .addoption(ECFOption().setname("STRUMPACK").setdescription("STRUMPACK"))
        .addoption(ECFOption().setname("PASTIX").setdescription("PaStiX"))
        .addoption(ECFOption().setname("SUPERLU_DIST").setdescription("SuperLU_dist"))
    );
}
