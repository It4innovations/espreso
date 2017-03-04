
#ifndef SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_ELASTICITY3D_ASSEMBLER_H_
#define SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_ELASTICITY3D_ASSEMBLER_H_

#include "../../../old_physics/linear/assembler.h"

namespace espreso {

struct LinearElasticity3DConfiguration;

struct LinearElasticity3D: public LinearPhysics
{
	LinearElasticity3D(Mesh &mesh, Constraints &constraints, const LinearElasticity3DConfiguration &configuration);

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1();
	void assembleB0();

	void saveMeshProperties(output::ResultStore &store);
	void saveMeshResults(output::ResultStore &store, const std::vector<std::vector<double> > &results);

	const LinearElasticity3DConfiguration &_configuration;

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_OLD_PHYSICS_LINEAR_ELASTICITY3D_ASSEMBLER_H_ */
