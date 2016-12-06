
#ifndef SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct PrecomputedPhysics: public Physics {

	virtual bool singular() const
	{
		return true;
	}

	virtual void assembleStiffnessMatrices()
	{
		ESINFO(PROGRESS2) << "Assemble matrices K and RHS.";
		#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);
			K[p].mtype = mtype;
			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	void saveMeshProperties(output::Store &store) { ESINFO(GLOBAL_ERROR) << "It is not possible to save mesh through API"; }
	void saveMeshResults(output::Store &store, const std::vector<std::vector<double> > &results) { ESINFO(GLOBAL_ERROR) << "It is not possible to save results through API"; }

	PrecomputedPhysics(
			APIMesh &mesh,
			Constraints &constraints,
			SparseMatrix::MatrixType mtype,
			const std::vector<Property> elementDOFs,
			const std::vector<Property> faceDOFs,
			const std::vector<Property> edgeDOFs,
			const std::vector<Property> pointDOFs,
			const std::vector<Property> midPointDOFs,
			double *rhs, size_t rhs_size)
	: Physics(mesh, constraints, mtype, elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs), _apimesh(mesh), _rhs(rhs), _rhs_size(rhs_size) {};
	virtual ~PrecomputedPhysics() {};

protected:
	virtual void composeSubdomain(size_t subdomain) =0;

	APIMesh &_apimesh;
	double *_rhs;
	size_t _rhs_size;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_ */
