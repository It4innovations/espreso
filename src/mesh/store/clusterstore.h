
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

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);
};

}

#endif /* SRC_MESH_STORE_CLUSTERSTORE_H_ */
