
#ifndef SOLVER_SPARSE_CPU_CLUSTERCPU_H_
#define SOLVER_SPARSE_CPU_CLUSTERCPU_H_

#include "../cluster.h"

class ClusterCPU: public ClusterBase
{

public:
	// Constructor
	ClusterCPU(eslocal cluster_index): ClusterBase(cluster_index) { };
	ClusterCPU(): ClusterBase() {};

	void Create_SC_perDomain( bool USE_FLOAT );
};


#endif /* SOLVER_SPARSE_CPU_CLUSTERCPU_H_ */
