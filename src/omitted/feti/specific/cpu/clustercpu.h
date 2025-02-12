
#ifndef SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_
#define SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_

#include "feti/specific/cluster.h"

namespace espreso {

class ClusterCPU: public ClusterBase
{

public:
    // Constructor
//    ClusterCPU(const FETISolverConfiguration &configuration, esint cluster_index): ClusterBase(configuration, cluster_index) { };
//    ClusterCPU(const FETISolverConfiguration &configuration): ClusterBase(configuration) {};
    ClusterCPU(const FETIConfiguration &configuration, DataHolder *instance_in): ClusterBase(configuration, instance_in) {};

    void Create_SC_perDomain( bool USE_FLOAT );
    void Create_Kinv_perDomain();
    void CreateDirichletPrec( DataHolder *instance );
    void SetupKsolvers ( );
};

}


#endif /* SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_ */
