
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct LinearElasticity: public LinearPhysics
{
	LinearElasticity(const Mesh &mesh)
	: LinearPhysics(mesh, { DOFType::DISPLACEMENT_X, DOFType::DISPLACEMENT_Y, DOFType::DISPLACEMENT_Z }, SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY_ASSEMBLER_H_ */
