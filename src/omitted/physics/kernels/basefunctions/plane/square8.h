
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_SQUARE8_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_SQUARE8_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Square8: public Element {

	static void setBaseFunctions(Element &self);
	void setGaussPointsForOrder(int order);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_SQUARE8_H_ */
