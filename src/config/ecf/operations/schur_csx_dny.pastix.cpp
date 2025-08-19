
#include "config/ecf/operations/schur_csx_dny.pastix.h"

#include "config/configuration.hpp"

using namespace espreso;

SchurCsxDnyPastixConfig::SchurCsxDnyPastixConfig()
{
    use_gpu = AUTOBOOL::AUTO;
    REGISTER(use_gpu, ECFMetaData()
        .setdescription({ "schur_csx_dny.pastix operation will internally use GPU" })
        .setdatatype({ ECFDataType::OPTION })
        .addoption(ECFOption().setname("AUTO").setdescription("Automatic selection"))
        .addoption(ECFOption().setname("TRUE").setdescription("GPU is enabled"))
        .addoption(ECFOption().setname("FALSE").setdescription("GPU is disabled"))
    );
}
