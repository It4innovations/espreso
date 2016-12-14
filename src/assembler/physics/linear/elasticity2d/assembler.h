
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_

#include "../assembler.h"
#include "../../../../config/linearelasticity2d.h"

namespace espreso {

struct LinearElasticity2D: public LinearPhysics
{
	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	LinearElasticity2D(Mesh &mesh, Constraints &constraints, const LinearElasticity2DConfiguration &configuration)
	: LinearPhysics(
			mesh, constraints, configuration.espreso,
			SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs),
	  _configuration(configuration) {};

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1();
	void assembleB0();

	void saveMeshProperties(store::Store &store);
	void saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results);

	static ELEMENT_BEHAVIOUR elementBehaviour;
	static Point angularVelocity;
	const LinearElasticity2DConfiguration &_configuration;

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

protected:
	void composeSubdomain(size_t subdomain);
};

}





#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_ */
