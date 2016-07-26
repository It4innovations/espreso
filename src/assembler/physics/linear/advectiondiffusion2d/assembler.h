
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion2D: public LinearPhysics
{
	AdvectionDiffusion2D(const Mesh &mesh)
	: LinearPhysics(mesh, { Property::TEMPERATURE }, { Property::FIXED_TEMPERATURE }, SparseMatrix::MatrixType::REAL_UNSYMMETRIC) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_ */
