
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION2D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion2D: public LinearPhysics
{
	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	AdvectionDiffusion2D(Mesh &mesh, Constraints &constraints, const ESPRESOSolver &configuration)
	: LinearPhysics(
			mesh, constraints, configuration,
			SparseMatrix::MatrixType::REAL_UNSYMMETRIC,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs) {};

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1();
	void assembleB0();

	void saveMeshProperties(store::Store &store);
	void saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results);

	static double sigma;
	static STABILIZATION stabilization;

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
