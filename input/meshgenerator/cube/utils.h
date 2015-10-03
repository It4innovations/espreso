
#ifndef INPUT_MESHGENERATOR_CUBE_UTILS_H_
#define INPUT_MESHGENERATOR_CUBE_UTILS_H_

#include "../../meshgenerator/cube/settings.h"

namespace esinput {

template <class TElement>
class Utils {

public:
	static void globalNodesCount(const CubeSettings &settings, esglobal nodes[]);
	static void clusterNodesCount(const CubeSettings &settings, eslocal nodes[]);
	static eslocal clusterElementsCount(const CubeSettings &settings);
};

}


#include "../../meshgenerator/cube/utils.hpp"


#endif /* INPUT_MESHGENERATOR_CUBE_UTILS_H_ */
