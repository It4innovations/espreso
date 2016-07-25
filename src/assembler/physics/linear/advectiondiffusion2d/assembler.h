
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion2D: public LinearPhysics
{
	AdvectionDiffusion2D(const Mesh &mesh)
	: LinearPhysics(mesh, { DOFType::TEMPERATURE }, SparseMatrix::MatrixType::REAL_UNSYMMETRIC)
	{
		if (mesh.initialConditions(InitialConditionType::HEAT_SOURCE) == NULL) {
			ESINFO(GLOBAL_ERROR) << "Set 'HEAT_SOURCES'";
		}
		if (mesh.initialConditions(InitialConditionType::TRANSLATION_MOTION_X) == NULL) {
			ESINFO(GLOBAL_ERROR) << "Set 'TRANSITION_MOTIONS' parameter 'x'";
		}
		if (mesh.initialConditions(InitialConditionType::TRANSLATION_MOTION_Y) == NULL) {
			ESINFO(GLOBAL_ERROR) << "Set 'TRANSITION_MOTIONS' parameter 'x'";
		}
	};

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_ */
