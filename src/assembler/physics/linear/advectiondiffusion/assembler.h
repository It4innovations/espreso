
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion: public LinearPhysics
{
	AdvectionDiffusion(const Mesh &mesh): LinearPhysics(mesh, 1) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION_ASSEMBLER_H_ */
