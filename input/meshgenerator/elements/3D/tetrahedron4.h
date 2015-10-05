
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_

#include "esmesh.h"

namespace esinput {

class Tetrahedron4 {

public:
	static size_t subelements;
	static size_t subnodes[3];

	static void addElements(std::vector<mesh::Element*> &elements, const eslocal indices[]);

};

}



#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_ */
