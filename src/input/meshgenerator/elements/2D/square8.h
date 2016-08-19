
#ifndef INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_
#define INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_

#include "esmesh.h"

namespace espreso {
namespace input {

class Square8 {

public:
	static size_t subelements;
	static size_t subnodes[3];

	static void addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[]);
	static void addFaces(std::vector<Element*> &faces, const eslocal indices[], size_t face) { ESINFO(GLOBAL_ERROR) << "Generator: plane element has no faces."; }

};

}
}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_ */
