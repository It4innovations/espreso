/*
 * clusterGPU.h
 *
 *  Created on: Feb 24, 2016
 *      Author: lriha
 */

#ifndef SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_


#include "../cluster.h"

namespace espreso {

class ClusterGPU: public ClusterBase
{

public:
	// Constructor
	ClusterGPU(eslocal cluster_index): ClusterBase(cluster_index) {};
	ClusterGPU(): ClusterBase() {};

	void Create_SC_perDomain( bool USE_FLOAT );
	void GetSchurComplement( bool USE_FLOAT, eslocal i );

	void SetupKsolvers ( );

	void multKplusGlobal_GPU   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );
};

}



#endif /* SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_ */
