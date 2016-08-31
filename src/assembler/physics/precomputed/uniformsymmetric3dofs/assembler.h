
#ifndef SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_UNIFORMSYMMETRIC3DOFS_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_UNIFORMSYMMETRIC3DOFS_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct UniformSymmetric3DOFs: public PrecomputedPhysics
{
	UniformSymmetric3DOFs(APIMesh &mesh, Constraints &constraints, double *rhs, size_t rhs_size)
	: PrecomputedPhysics(
			mesh, constraints,
			SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs,
			rhs, rhs_size) {};

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe);
	void assembleGluingMatrices() {};

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_UNIFORMSYMMETRIC3DOFS_ASSEMBLER_H_ */
