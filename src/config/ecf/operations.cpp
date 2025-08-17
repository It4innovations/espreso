
#include "config/ecf/operations.h"

#include "config/configuration.hpp"

using namespace espreso;

Operations::Operations()
{
    REGISTER(schur_csx_dny_spsolver, ECFMetaData()
            .setdescription({ "schur_csx_dny_spsolver config" }));

    REGISTER(schur_csx_dny_tria, ECFMetaData()
            .setdescription({ "schur_csx_dny_tria config" }));

    REGISTER(trsm_csx_dny_tria, ECFMetaData()
            .setdescription({ "trsm_csx_dny_tria config" }));

    REGISTER(herk_dnx_dny_tria, ECFMetaData()
            .setdescription({ "herk_dnx_dny_tria config" }));

    REGISTER(gpu_schur_hcsx_ddny_tria, ECFMetaData()
            .setdescription({ "gpu_schur_hcsx_ddny_tria config" }));

    REGISTER(gpu_trsm_hcsx_ddny_tria, ECFMetaData()
            .setdescription({ "gpu_trsm_hcsx_ddny_tria config" }));

    REGISTER(gpu_herk_ddnx_ddny_tria, ECFMetaData()
            .setdescription({ "gpu_herk_ddnx_ddny_tria config" }));
}
