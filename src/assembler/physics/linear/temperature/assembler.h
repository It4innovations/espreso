
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct Temperature: public LinearPhysics
{
	Temperature(const Mesh &mesh): LinearPhysics(mesh, 1) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_ */
