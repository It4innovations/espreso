
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_UTILS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_UTILS_H_

#include "settings.h"
#include "../utils.h"

namespace esinput {

template <class TElement>
class SphereUtils {

public:
	static eslocal surfaceNodesCount(const SphereSettings &settings);
	static eslocal ringNodesCount(const SphereSettings &settings);
};

}

#include "utils.hpp"


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_UTILS_H_ */
