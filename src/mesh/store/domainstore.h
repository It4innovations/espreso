
#ifndef SRC_MESH_STORE_DOMAINSTORE_H_
#define SRC_MESH_STORE_DOMAINSTORE_H_

#include "info.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct DomainStore: DistributedDataInfo {

	std::vector<size_t> distribution;

	std::vector<esint> gatherProcDistribution();

	serializededata<esint, esint>* nodes;
	std::vector<esint> elements;
	std::vector<int> cluster;

	serializededata<esint, esint> *dual;
	serializededata<esint, esint> *localDual;

	DomainStore();
	~DomainStore();
};

}

#endif /* SRC_MESH_STORE_DOMAINSTORE_H_ */
