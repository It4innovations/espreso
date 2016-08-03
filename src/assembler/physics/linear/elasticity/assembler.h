
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct LinearElasticity: public LinearPhysics
{
	LinearElasticity(const Mesh &mesh): LinearPhysics(mesh, 3) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY_ASSEMBLER_H_ */
