
#ifndef INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_
#define INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_

#include "quadraticplane.h"

namespace espreso {

struct Square8Generator: public QuadraticPlaneGenerator {

    Square8Generator();

    void pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const;
};

}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_ */
