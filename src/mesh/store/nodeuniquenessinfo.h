
#ifndef SRC_MESH_STORE_NODEUNIQUENESSINFO_H_
#define SRC_MESH_STORE_NODEUNIQUENESSINFO_H_

#include <vector>

namespace espreso {

struct NodeUniquenessInfo {
	esint nhalo;
	esint offset;
	esint size;
	esint totalSize;

	std::vector<esint> position;

	NodeUniquenessInfo(): nhalo(0), offset(0), size(0), totalSize(0) {}
};

}

#endif /* SRC_MESH_STORE_NODEUNIQUENESSINFO_H_ */
