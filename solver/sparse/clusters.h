
#ifndef SOLVER_SPARSE_CLUSTERS_H_
#define SOLVER_SPARSE_CLUSTERS_H_

#if defined(SOLVER_MKL)
#include "cpu/clustercpu.h"
	typedef ClusterCPU	Cluster;


#elif defined(SOLVER_PARDISO)
#include "cpu/clustercpu.h"
	typedef ClusterCPU	Cluster;

#elif defined(SOLVER_MUMPS)
#include "cpu/clustercpu.h"
	typedef ClusterCPU	Cluster;

#elif defined(SOLVER_MIC)
#include "acc/clusteracc.h"
	typedef ClusterAcc	Cluster;

#elif defined(SOLVER_CUDA)
#include "acc/clusteracc.h"
	typedef ClusterAcc	Cluster;


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif



#endif /* SOLVER_SPARSE_CLUSTERS_H_ */
