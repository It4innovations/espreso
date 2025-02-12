
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE3_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE3_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Triangle3: public Element {

    static void setBaseFunctions(Element &self);
    void setGaussPointsForOrder(int order);

    static void computeReferenceCoords(const MatrixDense & vertices, const MatrixDense & points, MatrixDense & result);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE3_H_ */
