
#ifndef INPUT_MESHGENERATOR_SPHERE_UTILS_H_
#define INPUT_MESHGENERATOR_SPHERE_UTILS_H_

#include "../../meshgenerator/sphere/settings.h"

namespace esinput {

template <class TElement>
class Utils {

public:
	static void globalNodesCount(const SphereSettings &settings, esglobal nodes[]);
	static void clusterNodesCount(const SphereSettings &settings, eslocal nodes[]);
	static eslocal clusterElementsCount(const SphereSettings &settings);
};

}


#include "../../meshgenerator/sphere/utils.hpp"


#endif /* INPUT_MESHGENERATOR_SPHERE_UTILS_H_ */
