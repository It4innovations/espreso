
#ifndef SRC_MESH_STORE_ELEMENTINFO_H_
#define SRC_MESH_STORE_ELEMENTINFO_H_

#include "info.h"
#include "mesh/element.h"

#include <vector>

namespace espreso {

struct ElementsDistributionInfo {

	std::vector<size_t> threads;
	DistributedDataInfo process;
	std::vector<DistributionInfo> code;

	ElementsDistributionInfo(): threads({ 0, 0 }), code(static_cast<int>(Element::CODE::SIZE)) {}
};

}

#endif /* SRC_MESH_STORE_ELEMENTINFO_H_ */
