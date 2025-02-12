
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_POINT_POINT1_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_POINT_POINT1_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Point1: public Element {

    static void setBaseFunctions(Element &self);
    void setGaussPointsForOrder(int order);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_POINT_POINT1_H_ */
