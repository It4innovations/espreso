
#ifndef INPUT_MESHGENERATOR_SPHERE_UTILS_H_
#define INPUT_MESHGENERATOR_SPHERE_UTILS_H_

#include "settings.h"

namespace esinput {

template <class TElement>
class SphereUtils {

public:
	static void clusterNodesCount(const SphereSettings &settings, eslocal nodes[]);
	static eslocal clusterElementsCount(const SphereSettings &settings);
};

}


#include "utils.hpp"


#endif /* INPUT_MESHGENERATOR_SPHERE_UTILS_H_ */
