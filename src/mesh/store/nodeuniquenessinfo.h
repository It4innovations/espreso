
#ifndef SRC_MESH_STORE_NODEUNIQUENESSINFO_H_
#define SRC_MESH_STORE_NODEUNIQUENESSINFO_H_

#include "basis/containers/point.h"

#include <vector>
#include <limits>

namespace espreso {

struct NodeUniquenessInfo {
	esint nhalo;
	esint offset;
	esint size;
	esint totalSize;

	Point min, max;

	std::vector<esint> position;

	NodeUniquenessInfo(): nhalo(0), offset(0), size(0), totalSize(0), min(std::numeric_limits<double>::max()), max(-min) {}
};

}

#endif /* SRC_MESH_STORE_NODEUNIQUENESSINFO_H_ */
