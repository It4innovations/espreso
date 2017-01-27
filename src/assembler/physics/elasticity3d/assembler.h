
#ifndef SRC_ASSEMBLER_PHYSICS_ELASTICITY3D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_ELASTICITY3D_ASSEMBLER_H_

#include "../assembler.h"
#include "../../../config/linearelasticity3d.h"

namespace espreso {

struct Elasticity3D: public Physics
{
	Elasticity3D(Mesh &mesh, Constraints &constraints, const LinearElasticity3DConfiguration &configuration)
	: Physics(
			mesh, constraints, configuration.espreso,
			MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs),
	  _configuration(configuration) {};

	bool singular() const
	{
		return true;
	}

	void assembleStiffnessMatrices();

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1();
	void assembleB0();

	void postProcess(store::Store &store, const std::vector<std::vector<double> > &solution);

	void saveMeshProperties(store::Store &store);
	void saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results);

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


#endif /* SRC_ASSEMBLER_PHYSICS_ELASTICITY3D_ASSEMBLER_H_ */
