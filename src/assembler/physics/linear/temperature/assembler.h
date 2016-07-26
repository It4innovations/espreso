
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct Temperature: public LinearPhysics
{
	Temperature(const Mesh &mesh)
	: LinearPhysics(mesh, { Property::TEMPERATURE }, { Property::FIXED_TEMPERATURE }, SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_ */
