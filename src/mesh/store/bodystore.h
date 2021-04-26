
#ifndef SRC_MESH_STORE_BODYSTORE_H_
#define SRC_MESH_STORE_BODYSTORE_H_

#include "info.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct BodyStore: UniqueDataInfo {

	BodyStore();
	~BodyStore();
};

}

#endif /* SRC_MESH_STORE_BODYSTORE_H_ */
