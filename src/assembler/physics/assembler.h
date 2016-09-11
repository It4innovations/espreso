
#ifndef SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_

#include "esmesh.h"
#include "essolver.h"
#include "esoutput.h"
#include "../constraints/constraints.h"

namespace espreso {

struct Physics {

	virtual bool singular() const =0;

	virtual void prepareMeshStructures() =0;
	virtual void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe) =0;
	virtual void assembleStiffnessMatrices() =0;
	virtual void assembleGluingMatrices() =0;

	virtual void saveMeshProperties(output::Store &store) =0;
	virtual void saveMeshResults(output::Store &store, const std::vector<std::vector<double> > &results) =0;
	virtual void saveMatrices()
	{
		ESINFO(PROGRESS2) << "Save matrices K and RHS.";
		for (size_t p = 0; p < K.size(); p++) {
			std::ofstream osK(Logging::prepareFile(p, "Kreg").c_str());
			osK << K[p];
			osK.close();
		}

		for (size_t p = 0; p < f.size(); p++) {
			std::ofstream osF(Logging::prepareFile(p, "f").c_str());
			osF << f[p];
			osF.close();
		}

		for (size_t p = 0; p < R1.size(); p++) {
			std::ofstream osR(Logging::prepareFile(p, "R1_").c_str());
			SparseMatrix tmpR = R1[p];
			tmpR.ConvertDenseToCSR(0);
			osR << tmpR;
			osR.close();
		}

		for (size_t p = 0; p < R2.size(); p++) {
			std::ofstream osR(Logging::prepareFile(p, "R2_").c_str());
			SparseMatrix tmpR = R2[p];
			tmpR.ConvertDenseToCSR(0);
			osR << tmpR;
			osR.close();
		}
	}

	std::vector<Property> elementDOFs;
	std::vector<Property> faceDOFs;
	std::vector<Property> edgeDOFs;
	std::vector<Property> pointDOFs;
	std::vector<Property> midPointDOFs;

	std::vector<size_t> matrixSize;

	SparseMatrix::MatrixType mtype;
	std::vector<SparseMatrix> K, R1, R2, RegMat;
	std::vector<std::vector<double> > f;

	Physics(Mesh &mesh,
			Constraints &constraints,
			SparseMatrix::MatrixType mtype,
			const std::vector<Property> elementDOFs,
			const std::vector<Property> faceDOFs,
			const std::vector<Property> edgeDOFs,
			const std::vector<Property> pointDOFs,
			const std::vector<Property> midPointDOFs):
	_mesh(mesh),
	_constraints(constraints),
	mtype(mtype),
	elementDOFs(elementDOFs),
	faceDOFs(faceDOFs),
	edgeDOFs(edgeDOFs),
	pointDOFs(pointDOFs),
	midPointDOFs(midPointDOFs),
	matrixSize(mesh.parts()) {};

	virtual ~Physics() {};

protected:
	Mesh& _mesh;
	Constraints &_constraints;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_ */
