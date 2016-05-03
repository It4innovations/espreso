
#ifndef SOLVER_SPECIFIC_CLUSTERS_H_
#define SOLVER_SPECIFIC_CLUSTERS_H_

#if defined(SOLVER_MKL)
#include "cpu/clustercpu.h"

namespace espreso {
	typedef ClusterCPU	Cluster;
}


#elif defined(SOLVER_PARDISO)
#include "cpu/clustercpu.h"

namespace espreso {
	typedef ClusterCPU	Cluster;
}

#elif defined(SOLVER_MUMPS)
#include "cpu/clustercpu.h"

namespace espreso {
	typedef ClusterCPU	Cluster;
}

#elif defined(SOLVER_MIC)
#include "acc/clusteracc.h"

namespace espreso {
	typedef ClusterAcc	Cluster;
}

#elif defined(SOLVER_CUDA)
#include "acc/clusterGPU.h"

namespace espreso {
	typedef ClusterGPU	Cluster;
}

#elif defined(SOLVER_CUDA_7)
#include "cpu/clustercpu.h"

namespace espreso {
	typedef ClusterCPU	Cluster;
}


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif



#endif /* SOLVER_SPECIFIC_CLUSTERS_H_ */
