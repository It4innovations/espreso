
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion2D: public LinearPhysics
{
	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	AdvectionDiffusion2D(const Mesh &mesh)
	: LinearPhysics(mesh, { Property::TEMPERATURE }, SparseMatrix::MatrixType::REAL_UNSYMMETRIC) {};

	static double sigma;
	static STABILIZATION stabilization;


protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_ */
