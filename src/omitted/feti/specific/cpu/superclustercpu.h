
#ifndef SRC_SOLVER_SPECIFIC_CPU_SUPERCLUSTERCPU_H_
#define SRC_SOLVER_SPECIFIC_CPU_SUPERCLUSTERCPU_H_

#include "feti/specific/supercluster.h"

namespace espreso {

class SuperClusterCPU : public SuperClusterBase
{
    public:

    SuperClusterCPU( const FETIConfiguration & configuration, DataHolder *instance_in ):
        SuperClusterBase( configuration, instance_in ) { }

    void Create_SC_perDomain(bool USE_FLOAT) {
        //bool
        USE_FLOAT = false;
        if (configuration.schur_precision == FETIConfiguration::FLOAT_PRECISION::SINGLE) {
            USE_FLOAT = true;
        }

        for (size_t c = 0; c < clusters.size(); c++) {
            clusters[c].Create_SC_perDomain(USE_FLOAT);
        }
    }

    void SetupKsolvers () {
        for (size_t c = 0; c < clusters.size(); c++) {
            clusters[c].SetupKsolvers();
        }
    }
    
    void SetupPreconditioner() {
        for (size_t c = 0; c < clusters.size(); c++) {
            clusters[c].SetupPreconditioner();
        }
    }

};

}

#endif
