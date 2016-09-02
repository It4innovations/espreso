
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_

#include "esmesh.h"
#include "../element.h"

namespace espreso {
namespace input {

class Pyramid5 {

public:
	static size_t subelements;
	static size_t subnodes[3];

	static void addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[]);
	static void addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFaces face);
	static void addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdges edge) { ESINFO(GLOBAL_ERROR) << "generator Edge"; };
	static void pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeEdges edge) { ESINFO(GLOBAL_ERROR) << "generator Nodes"; };
	static void pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFaces face);

};

}
}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_ */
