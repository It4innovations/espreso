
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_TETRAHEDRON4_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_TETRAHEDRON4_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Tetrahedron4: public Element {

    static void setBaseFunctions(Element &self);
    void setGaussPointsForOrder(int order);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_TETRAHEDRON4_H_ */
