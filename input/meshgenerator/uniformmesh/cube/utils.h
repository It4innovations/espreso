
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_UTILS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_UTILS_H_

#include "settings.h"
#include "../utils.h"

namespace esinput {

template <class TElement>
class CubeUtils {

public:
	static void globalNodesCount(const CubeSettings &settings, esglobal nodes[]);
};

}

#include "utils.hpp"


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_UTILS_H_ */
