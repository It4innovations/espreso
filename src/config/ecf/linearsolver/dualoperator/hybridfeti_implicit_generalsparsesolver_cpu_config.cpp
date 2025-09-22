
#include "config/ecf/linearsolver/dualoperator/hybridfeti_implicit_generalsparsesolver_cpu_config.h"

#include "config/configuration.hpp"

using namespace espreso;

DualopHybridfetiImplicitGeneralSparseSolverCpuConfig::DualopHybridfetiImplicitGeneralSparseSolverCpuConfig()
{
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
