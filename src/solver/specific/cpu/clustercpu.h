
#ifndef SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_
#define SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_

#include "../cluster.h"

namespace espreso {

class ClusterCPU: public ClusterBase
{

public:
	// Constructor
	ClusterCPU(eslocal cluster_index): ClusterBase(cluster_index) { };
	ClusterCPU(): ClusterBase() {};

	void Create_SC_perDomain( bool USE_FLOAT );
    void Create_Kinv_perDomain();
	void SetupKsolvers ( );
};

}


#endif /* SOLVER_SPECIFIC_CPU_CLUSTERCPU_H_ */
