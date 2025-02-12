
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_HEXAHEDRON20_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_HEXAHEDRON20_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Hexahedron20: public Element {

    static void setBaseFunctions(Element &self);
    void setGaussPointsForOrder(int order);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_HEXAHEDRON20_H_ */
