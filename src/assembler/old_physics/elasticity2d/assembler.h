
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct LinearElasticity2DConfiguration;

struct Elasticity2D: public Physics
{
	Elasticity2D(Mesh &mesh, Constraints &constraints, const LinearElasticity2DConfiguration &configuration);

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

	void saveMeshProperties(store::ResultStore &store);
	void saveMeshResults(store::ResultStore &store, const std::vector<std::vector<double> > &results);

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
