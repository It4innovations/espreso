
#ifndef SRC_MESH_STORE_DOMAINSTORE_H_
#define SRC_MESH_STORE_DOMAINSTORE_H_

#include "info.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct DomainStore: UniqueDataInfo {

	std::vector<size_t> distribution;

	std::vector<esint> gatherProcDistribution();

	DomainStore();
	~DomainStore();
};

}

#endif /* SRC_MESH_STORE_DOMAINSTORE_H_ */
