
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_UTILS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_UTILS_H_

#include "settings.h"

namespace espreso {
namespace input {

template <class TElement>
class UniformUtils {

public:
	static void clusterNodesCount(const UniformSettings &settings, eslocal nodes[]);
	static eslocal clusterElementsCount(const UniformSettings &settings);
};

}
}


#include "utils.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_UTILS_H_ */
