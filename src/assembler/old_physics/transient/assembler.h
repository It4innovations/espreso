
#ifndef SRC_ASSEMBLER_OLD_PHYSICS_TRANSIENT_ASSEMBLER_H_
#define SRC_ASSEMBLER_OLD_PHYSICS_TRANSIENT_ASSEMBLER_H_

#include "../../old_physics/assembler.h"

namespace espreso {

struct TransientPhysics: public Physics {

	virtual bool singular() const
	{
		return false;
	}

	virtual void assembleStiffnessMatrices();

	virtual void saveStiffnessMatrices();

	TransientPhysics(
			Mesh &mesh,
			Constraints &constraints,
			const ESPRESOSolver &configuration,
			MatrixType mtype,
			const std::vector<Property> elementDOFs,
			const std::vector<Property> faceDOFs,
			const std::vector<Property> edgeDOFs,
			const std::vector<Property> pointDOFs,
			const std::vector<Property> midPointDOFs);

	virtual ~TransientPhysics() {}

	std::vector<SparseMatrix> M;
	std::vector<double> A;

protected:
	virtual void composeSubdomain(size_t subdomain) =0;
};

}


#endif /* SRC_ASSEMBLER_OLD_PHYSICS_TRANSIENT_ASSEMBLER_H_ */
