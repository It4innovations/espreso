
#ifndef INPUT_ANSYS_UTILS_H_
#define INPUT_ANSYS_UTILS_H_

#include <cstdlib>

namespace espreso {

class Element;

namespace input {

class AnsysUtils {
public:
	static Element* createElement(eslocal *indices, eslocal n, eslocal *params);
	static Element* createElement(eslocal *indices, eslocal n, eslocal *params, int eType);
};


}
}




#endif /* INPUT_ANSYS_UTILS_H_ */
