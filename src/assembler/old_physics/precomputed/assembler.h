
#ifndef SRC_ASSEMBLER_OLD_PHYSICS_PRECOMPUTED_ASSEMBLER_H_
#define SRC_ASSEMBLER_OLD_PHYSICS_PRECOMPUTED_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

class APIMesh;

struct PrecomputedPhysics: public OldPhysics {

	virtual bool singular() const
	{
		return true;
	}

	virtual void assembleStiffnessMatrices();

	void saveMeshProperties(output::ResultStore &store);
	void saveMeshResults(output::ResultStore &store, const std::vector<std::vector<double> > &results);

	PrecomputedPhysics(
			APIMesh &mesh,
			Constraints &constraints,
			const ESPRESOSolver &configuration,
			MatrixType mtype,
			const std::vector<Property> elementDOFs,
			const std::vector<Property> faceDOFs,
			const std::vector<Property> edgeDOFs,
			const std::vector<Property> pointDOFs,
			const std::vector<Property> midPointDOFs,
			double *rhs, size_t rhs_size);
	virtual ~PrecomputedPhysics() {};

protected:
	virtual void composeSubdomain(size_t subdomain) =0;

	APIMesh &_apimesh;
	double *_rhs;
	size_t _rhs_size;
};

}


#endif /* SRC_ASSEMBLER_OLD_PHYSICS_PRECOMPUTED_ASSEMBLER_H_ */
