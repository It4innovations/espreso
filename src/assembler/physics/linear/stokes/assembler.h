
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_STOKES_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_STOKES_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct Stokes: public LinearPhysics
{
	Stokes(const Mesh &mesh): LinearPhysics(mesh, 3, SparseMatrix::MatrixType::REAL_SYMMETRIC_INDEFINITE) {};

protected:
	void composeSubdomain(size_t subdomain);
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_STOKES_ASSEMBLER_H_ */
