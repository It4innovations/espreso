
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA15_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA15_H_

#include "esmesh.h"

namespace espreso {
namespace input {

class Prisma15 {

public:
	static size_t subelements;
	static size_t subnodes[3];

	static void addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[]);

};

}
}



#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA15_H_ */
