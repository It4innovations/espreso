/*
 * clusterGPU.h
 *
 *  Created on: Feb 24, 2016
 *      Author: lriha
 */

#ifndef SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_


#include "../cluster.h"

class ClusterGPU: public ClusterBase
{

public:
	// Constructor
	ClusterGPU(eslocal cluster_index): ClusterBase(cluster_index) {};
	ClusterGPU(): ClusterBase() {};

	void Create_SC_perDomain( bool USE_FLOAT );
	void SetupKsolvers ( );
};



#endif /* SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_ */
