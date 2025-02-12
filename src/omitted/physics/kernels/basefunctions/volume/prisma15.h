
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_PRISMA15_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_PRISMA15_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Prisma15: public Element {

    static void setBaseFunctions(Element &self);
    void setGaussPointsForOrder(int order);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_PRISMA15_H_ */
