
#include "config/ecf/linearsolver/dirichlet_generalschur_config.h"

#include "config/configuration.hpp"

using namespace espreso;

DirichletGeneralSchurConfig::DirichletGeneralSchurConfig()
{
    parallel_set = AUTOBOOL::AUTO;
    REGISTER(parallel_set, ECFMetaData()
        .setdescription({ "Parallelism of the main loops in set" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Main loops will be parallel"))
        .addoption(ECFOption().setname("FALSE").setdescription("Main loops will be sequential"))
    );

    parallel_update = AUTOBOOL::AUTO;
    REGISTER(parallel_update, ECFMetaData()
        .setdescription({ "Parallelism of the main loops in update" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Main loops will be parallel"))
        .addoption(ECFOption().setname("FALSE").setdescription("Main loops will be sequential"))
    );

    parallel_apply = AUTOBOOL::AUTO;
    REGISTER(parallel_apply, ECFMetaData()
        .setdescription({ "Parallelism of the main loops in apply" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Main loops will be parallel"))
        .addoption(ECFOption().setname("FALSE").setdescription("Main loops will be sequential"))
    );

    timers_outer = AUTOBOOL::AUTO;
    REGISTER(timers_outer, ECFMetaData()
        .setdescription({ "Verbosity of outer timers" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Outer timers will be printed"))
        .addoption(ECFOption().setname("FALSE").setdescription("Outer timers will not be printed"))
    );

    timers_inner = AUTOBOOL::AUTO;
    REGISTER(timers_inner, ECFMetaData()
        .setdescription({ "Verbosity of inner timers" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Inner timers will be printed"))
        .addoption(ECFOption().setname("FALSE").setdescription("Inner timers will not be printed"))
    );

    print_config = AUTOBOOL::AUTO;
    REGISTER(print_config, ECFMetaData()
        .setdescription({ "Verbosity of info method regarding config" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Config parameters will be printed"))
        .addoption(ECFOption().setname("FALSE").setdescription("Config parameters will not be printed"))
    );

    order_sc = MATRIX_ORDER::AUTO;
    REGISTER(order_sc, ECFMetaData()
        .setdescription({ "Memory order of the dense Dirichlet preconditioner matrix" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    schur_impl_cpu = SCHUR_IMPL_CPU::AUTO;
    REGISTER(schur_impl_cpu, ECFMetaData()
        .setdescription({ "Schur wrapper to be used" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("MANUAL_SIMPLE").setdescription("Simple manual assembly"))
        .addoption(ECFOption().setname("TRIANGULAR").setdescription("Triangular"))
        .addoption(ECFOption().setname("MKLPARDISO").setdescription("MKL Pardiso"))
        .addoption(ECFOption().setname("SPARSE_SOLVER").setdescription("Manual using sparse solver wrapper"))
        .addoption(ECFOption().setname("MUMPS").setdescription("MUMPS"))
        .addoption(ECFOption().setname("PASTIX").setdescription("PaStiX"))
    );

    schur_impl_gpu = SCHUR_IMPL_GPU::AUTO;
    REGISTER(schur_impl_gpu, ECFMetaData()
        .setdescription({ "Schur wrapper to be used" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("MANUAL_SIMPLE").setdescription("Simple manual assembly"))
        .addoption(ECFOption().setname("TRIANGULAR").setdescription("Triangular"))
    );
}
