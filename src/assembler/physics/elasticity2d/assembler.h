
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ELASTICITY2D_ASSEMBLER_H_

#include "../assembler.h"
#include "../../../config/linearelasticity2d.h"

namespace espreso {

struct Elasticity2D: public Physics
{
	Elasticity2D(Mesh &mesh, Constraints &constraints, const LinearElasticity2DConfiguration &configuration)
	: Physics(
			mesh, constraints, configuration.espreso,
			SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
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

	void saveMeshProperties(store::Store &store);
	void saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results);

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
