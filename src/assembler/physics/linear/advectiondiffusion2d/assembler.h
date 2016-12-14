
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_

#include "../assembler.h"
#include "../../../../config/advectiondiffusion2d.h"

namespace espreso {

struct AdvectionDiffusion2D: public LinearPhysics
{
	AdvectionDiffusion2D(Mesh &mesh, Constraints &constraints, const AdvectionDiffusion2DConfiguration &configuration)
	: LinearPhysics(
			mesh, constraints, configuration.espreso,
			SparseMatrix::MatrixType::REAL_UNSYMMETRIC,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs),
	  _configuration(configuration) {};

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1();
	void assembleB0();

	void saveMeshProperties(store::Store &store);
	void saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results);

	const AdvectionDiffusion2DConfiguration &_configuration;

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_ */
