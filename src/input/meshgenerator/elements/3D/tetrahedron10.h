
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_

#include "esmesh.h"
#include "../element.h"

namespace espreso {
namespace input {

class Tetrahedron10 {

public:
	static size_t subelements;
	static size_t subnodes[3];

	static void addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[]);
	static void addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFaces face) { ESINFO(GLOBAL_ERROR) << "generator Face"; };
	static void addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdges edge) { ESINFO(GLOBAL_ERROR) << "generator Edge"; };

};

}
}



#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_ */
