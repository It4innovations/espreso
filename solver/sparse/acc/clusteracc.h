
#ifndef SOLVER_SPARSE_ACC_CLUSTERACC_H_
#define SOLVER_SPARSE_ACC_CLUSTERACC_H_

#include "../cluster.h"

class ClusterAcc: public ClusterBase
{

public:
	// Constructor
	ClusterAcc(eslocal cluster_index): ClusterBase(cluster_index) {};
	ClusterAcc(): ClusterBase() {};

	void Create_SC_perDomain( bool USE_FLOAT );
};



#endif /* SOLVER_SPARSE_ACC_CLUSTERACC_H_ */
