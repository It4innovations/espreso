
#ifndef SRC_MESH_STORE_CLUSTERSTORE_H_
#define SRC_MESH_STORE_CLUSTERSTORE_H_

#include "info.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct ClusterStore: DistributedDataInfo {

	ClusterStore();
	~ClusterStore();
};

}

#endif /* SRC_MESH_STORE_CLUSTERSTORE_H_ */
