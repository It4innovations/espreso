
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION3D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION3D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion3D: public LinearPhysics
{
	AdvectionDiffusion3D(const Mesh &mesh)
	: LinearPhysics(mesh, { Property::TEMPERATURE }, { Property::FIXED_TEMPERATURE }, SparseMatrix::MatrixType::REAL_UNSYMMETRIC) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION3D_ASSEMBLER_H_ */
