
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_SQUARE4_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_SQUARE4_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Square4: public Element {

	static void setBaseFunctions(Element &self);
	void setGaussPointsForOrder(int order);

	void computeReferenceCoords(const MatrixDense & vertices, const MatrixDense & points, MatrixDense & result);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_SQUARE4_H_ */
