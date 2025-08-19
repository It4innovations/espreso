
#include "config/ecf/operations/solver_csx.pastix.h"

#include "config/configuration.hpp"

using namespace espreso;

SolverCsxPastixConfig::SolverCsxPastixConfig()
{
    use_gpu = AUTOBOOL::AUTO;
    REGISTER(use_gpu, ECFMetaData()
        .setdescription({ "solver_csx.pastix operation will internally use GPU" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("GPU is enabled"))
        .addoption(ECFOption().setname("FALSE").setdescription("GPU is disabled"))
    );
}
