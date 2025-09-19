
#include "config/ecf/linearsolver/dualoperator/hybridfeti_explicit_generalschur_gpu_config.h"

#include "config/configuration.hpp"

using namespace espreso;

DualopHybridfetiExplicitGeneralSchurGpuConfig::DualopHybridfetiExplicitGeneralSchurGpuConfig()
{
    mainloop_update_split = MAINLOOP_UPDATE_SPLIT::AUTO;
    REGISTER(mainloop_update_split, ECFMetaData()
        .setdescription({ "Numerical factorization and Schur computation in update" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("COMBINED").setdescription("Combined in a single loop"))
        .addoption(ECFOption().setname("SEPARATE").setdescription("Separated into two loops"))
    );

    synchronize_after_update_mainloop = AUTOBOOL::AUTO;
    REGISTER(synchronize_after_update_mainloop, ECFMetaData()
        .setdescription({ "Synchronize with GPU after main loop completion" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Synchronize"))
        .addoption(ECFOption().setname("FALSE").setdescription("Don't synchronize"))
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

    order_F = MATRIX_ORDER::AUTO;
    REGISTER(order_F, ECFMetaData()
        .setdescription({ "Memory order of the dense matrix F" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major"))
        .addoption(ECFOption().setname("COL_MAJOR").setdescription("Column-major"))
    );

    schur_impl = SCHUR_IMPL::AUTO;
    REGISTER(schur_impl, ECFMetaData()
        .setdescription({ "Schur GPU wrapper to be used" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("MANUAL_SIMPLE").setdescription("Simple manual assembly"))
        .addoption(ECFOption().setname("TRIANGULAR").setdescription("Triangular"))
    );

    apply_where = CPU_GPU::AUTO;
    REGISTER(apply_where, ECFMetaData()
        .setdescription({ "Where does the explicit apply occur" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("CPU").setdescription("On the CPU"))
        .addoption(ECFOption().setname("GPU").setdescription("On the GPU"))
    );

    apply_wait_intermediate = AUTOBOOL::AUTO;
    REGISTER(apply_wait_intermediate, ECFMetaData()
        .setdescription({ "Wait for GPU in apply also before the external function" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("Yes, wait"))
        .addoption(ECFOption().setname("FALSE").setdescription("No, let GPU execute asynchronously"))
    );
}
