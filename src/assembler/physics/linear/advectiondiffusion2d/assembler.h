
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion2DConfiguration;

struct AdvectionDiffusion2D: public LinearPhysics
{
	AdvectionDiffusion2D(Mesh &mesh, Constraints &constraints, const AdvectionDiffusion2DConfiguration &configuration);

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1();
	void assembleB0();

	void postProcess(store::ResultStore &store, const std::vector<std::vector<double> > &solution);

	void saveMeshProperties(store::ResultStore &store);
	void saveMeshResults(store::ResultStore &store, const std::vector<std::vector<double> > &results);

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
