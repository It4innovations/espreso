
#ifndef SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_LAMESTEKLOVPOINCARE_ASSEMBLER_H_
#define SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_LAMESTEKLOVPOINCARE_ASSEMBLER_H_

#include "../../../old_physics/linear/assembler.h"

namespace espreso {

struct LameSteklovPoincare: public LinearPhysics
{
	LameSteklovPoincare(Mesh &mesh, Constraints &constraints, const ESPRESOSolver &configuration)
	: LinearPhysics(
			mesh, constraints, configuration,
			MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs) {};

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1() {};
	void assembleB0() {};

	void saveMeshProperties(store::ResultStore &store);
	void saveMeshResults(store::ResultStore &store, const std::vector<std::vector<double> > &results);

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_LAMESTEKLOVPOINCARE_ASSEMBLER_H_ */
