
#ifndef INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE4_H_
#define INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE4_H_

#include "esmesh.h"
#include "../element.h"

namespace espreso {
namespace input {

class Square4 {

public:
	static size_t subelements;
	static size_t subnodes[3];

	static void addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[]);
	static void addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFaces face) { ESINFO(GLOBAL_ERROR) << "Generator: plane element has no faces."; }
	static void addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdges edge);

};

}
}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE4_H_ */
