/*
 * superclusters.h
 *
 *  Created on: Jun 8, 2017
 *      Author: lriha
 */

#ifndef SRC_SOLVER_SPECIFIC_SUPERCLUSTERS_H_
#define SRC_SOLVER_SPECIFIC_SUPERCLUSTERS_H_

#if defined(SOLVER_MKL)
#include "cpu/superclustercpu.h"

namespace espreso {
    typedef SuperClusterCPU    SuperCluster;
}


#elif defined(SOLVER_PARDISO)
#include "cpu/superclustercpu.h"

namespace espreso {
    typedef SuperClusterCPU    SuperCluster;
}

#elif defined(SOLVER_MUMPS)
#include "cpu/superclustercpu.h"

namespace espreso {
    typedef SuperClusterCPU    SuperCluster;
}

#elif defined(SOLVER_MIC)
#include "acc/superclusteracc.h"

namespace espreso {
    typedef SuperClusterAcc    SuperCluster;
}

#elif defined(SOLVER_CUDA)
#include "cpu/superclustercpu.h"

namespace espreso {
    typedef SuperClusterCPU    SuperCluster;
}

#elif defined(SOLVER_CUDA_7)
#include "cpu/superclustercpu.h"

namespace espreso {
    typedef SuperClusterCPU    SuperCluster;
}

#elif defined(SOLVER_DISSECTION)
#include "cpu/superclustercpu.h"

namespace espreso {
    typedef SuperClusterCPU    SuperCluster;
}


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SRC_SOLVER_SPECIFIC_SUPERCLUSTERS_H_ */
