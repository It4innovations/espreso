
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE6_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE6_H_

#include "physics/kernels/basefunctions/basefunctions.h"

namespace espreso {

struct Triangle6: public Element {

	static void setBaseFunctions(Element &self);
	void setGaussPointsForOrder(int order);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE6_H_ */
