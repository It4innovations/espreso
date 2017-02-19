
#ifndef SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_
#define SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_

#include "../cluster.h"

namespace espreso {

class ClusterCPU: public ClusterBase
{

public:
	// Constructor
	ClusterCPU(const ESPRESOSolver &configuration, eslocal cluster_index): ClusterBase(configuration, cluster_index) { };
	ClusterCPU(const ESPRESOSolver &configuration): ClusterBase(configuration) {};

	void Create_SC_perDomain( bool USE_FLOAT );
    void Create_Kinv_perDomain();
    void CreateDirichletPrec( Instance *instance );
	void SetupKsolvers ( );
};

}


#endif /* SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_ */
