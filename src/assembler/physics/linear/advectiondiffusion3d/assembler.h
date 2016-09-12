
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION3D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION3D_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct AdvectionDiffusion3D: public LinearPhysics
{
	AdvectionDiffusion3D(Mesh &mesh, Constraints &constraints)
	: LinearPhysics(
			mesh, constraints,
			SparseMatrix::MatrixType::REAL_UNSYMMETRIC,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs) {};

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe);
	void makeStiffnessMatricesRegular();
	void assembleGluingMatrices() {};

	void saveMeshProperties(output::Store &store);
	void saveMeshResults(output::Store &store, const std::vector<std::vector<double> > &results);

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ADVECTIONDIFFUSION3D_ASSEMBLER_H_ */
