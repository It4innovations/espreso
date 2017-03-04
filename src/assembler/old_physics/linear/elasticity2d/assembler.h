
#ifndef SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_

#include "../../../old_physics/linear/assembler.h"

namespace espreso {

struct LinearElasticity2DConfiguration;

struct LinearElasticity2D: public LinearPhysics
{
	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	LinearElasticity2D(Mesh &mesh, Constraints &constraints, const LinearElasticity2DConfiguration &configuration);

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1();
	void assembleB0();

	void saveMeshProperties(output::ResultStore &store);
	void saveMeshResults(output::ResultStore &store, const std::vector<std::vector<double> > &results);

	static ELEMENT_BEHAVIOUR elementBehaviour;
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





#endif /* SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_ */
