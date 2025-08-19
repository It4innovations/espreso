
#ifndef SRC_CONFIG_ECF_OPERATIONS_H_
#define SRC_CONFIG_ECF_OPERATIONS_H_

#include "config/description.h"
#include "operations/schur_csx_dny.spsolver.h"
#include "operations/schur_csx_dny.tria.h"
#include "operations/trsm_csx_dny_tria.h"
#include "operations/herk_dnx_dny_tria.h"
#include "operations/gpu_schur_hcsx_ddny.tria.h"
#include "operations/gpu_trsm_hcsx_ddny_tria.h"
#include "operations/gpu_herk_ddnx_ddny_tria.h"
#include "operations/schur_csx_dny.pastix.h"
#include "operations/solver_csx.pastix.h"
#include "operations/solver_csx.superlu_dist.h"

namespace espreso {

struct Operations: public ECFDescription {

    SchurCsxDnySpsolverConfig schur_csx_dny_spsolver;
    SchurCsxDnyTriaConfig schur_csx_dny_tria;
    TrsmCsxDnyTriaConfig trsm_csx_dny_tria;
    HerkDnxDnyTriaConfig herk_dnx_dny_tria;

    GpuSchurHcsxDdnyTriaConfig gpu_schur_hcsx_ddny_tria;
    GpuTrsmHcsxDdnyTriaConfig gpu_trsm_hcsx_ddny_tria;
    GpuHerkDdnxDdnyTriaConfig gpu_herk_ddnx_ddny_tria;

    SchurCsxDnyPastixConfig schur_csx_dny_pastix;
    SolverCsxPastixConfig solver_csx_pastix;
    SolverCsxSuperludistConfig solver_csx_superlu_dist;

    Operations();

};

}

#endif /* SRC_CONFIG_ECF_OPERATIONS_H_ */
