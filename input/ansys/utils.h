
#ifndef INPUT_ANSYS_UTILS_H_
#define INPUT_ANSYS_UTILS_H_

#include <cstdlib>

#include "esmesh.h"

namespace esinput {

class AnsysUtils {
public:
	static mesh::Element* createElement(eslocal *indices, eslocal n);
};


}




#endif /* INPUT_ANSYS_UTILS_H_ */
