
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_

#include "esmesh.h"

namespace espreso {
namespace input {

class Tetrahedron4 {

public:
	static size_t subelements;
	static size_t subnodes[3];

	static void addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[]);
	static void addFaces(std::vector<Element*> &faces, const eslocal indices[], size_t face) { ESINFO(GLOBAL_ERROR) << "generator Face"; };

};

}
}



#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_ */
