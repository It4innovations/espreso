/*
 * superclusters.h
 *
 *  Created on: Jun 8, 2017
 *      Author: lriha
 */

#ifndef SRC_SOLVER_SPECIFIC_SUPERCLUSTERS_H_
#define SRC_SOLVER_SPECIFIC_SUPERCLUSTERS_H_

#if defined(SOLVER_MKL)
#include "supercluster.h"//"cpu/clustercpu.h"

namespace espreso {
	;//typedef ClusterCPU	Cluster;
}


#elif defined(SOLVER_PARDISO)
#include "supercluster.h"//"cpu/clustercpu.h"

namespace espreso {
	;//typedef ClusterCPU	Cluster;
}

#elif defined(SOLVER_MUMPS)
#include "supercluster.h"//"cpu/clustercpu.h"

namespace espreso {
	;//typedef ClusterCPU	Cluster;
}

#elif defined(SOLVER_MIC)
#include "supercluster.h"//"cpu/clustercpu.h"

namespace espreso {
	;//typedef ClusterCPU	Cluster;
}

#elif defined(SOLVER_CUDA)
#include "supercluster.h"//"cpu/clustercpu.h"

namespace espreso {
	;//typedef ClusterCPU	Cluster;
}

#elif defined(SOLVER_CUDA_7)
#include "supercluster.h"//"cpu/clustercpu.h"

namespace espreso {
	;//typedef ClusterCPU	Cluster;
}

#elif defined(SOLVER_DISSECTION)
#include "supercluster.h"//"cpu/clustercpu.h"

namespace espreso {
	;//typedef ClusterCPU	Cluster;
}


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SRC_SOLVER_SPECIFIC_SUPERCLUSTERS_H_ */
