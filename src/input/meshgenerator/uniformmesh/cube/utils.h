
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_UTILS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_UTILS_H_

#include "settings.h"
#include "../utils.h"

namespace espreso {
namespace input {

template <class TElement>
class CubeUtils {

public:
	static void globalNodesCount(const CubeSettings &settings, esglobal nodes[]);
	static void computeInterval(const CubeSettings &settings, size_t cluster[], const Interval &interval, size_t start[], size_t end[]);
	static CubeEdges cubeEdge(const CubeSettings &settings, size_t cluster[], const Interval &interval);
	static CubeFaces cubeFace(const CubeSettings &settings, size_t cluster[], const Interval &interval);
};

}
}

#include "utils.hpp"


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_UTILS_H_ */
