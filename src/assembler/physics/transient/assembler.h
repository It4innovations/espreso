
#ifndef SRC_ASSEMBLER_PHYSICS_TRANSIENT_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_TRANSIENT_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct TransientPhysics: public Physics {

	virtual bool singular() const
	{
		return false;
	}

	virtual void assembleStiffnessMatrices()
	{
		ESINFO(PROGRESS2) << "Assemble matrices K, M, and RHS.";
		#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);
			K[p].mtype = mtype;
			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	virtual void saveStiffnessMatrices()
	{
		ESINFO(PROGRESS2) << "Save matrices K, M, RHS, and A constant";
		for (size_t p = 0; p < _mesh.parts(); p++) {
			std::ofstream osK(Logging::prepareFile(p, "K").c_str());
			osK << K[p];
			osK.close();

			std::ofstream osM(Logging::prepareFile(p, "M").c_str());
			osM << M[p];
			osM.close();

			std::ofstream osF(Logging::prepareFile(p, "f").c_str());
			osF << f[p];
			osF.close();
		}

		std::ofstream osA(Logging::prepareFile("A").c_str());
		osA << A;
		osA.close();
	}

	TransientPhysics(
			Mesh &mesh,
			Constraints &constraints,
			const ESPRESOSolver &configuration,
			SparseMatrix::MatrixType mtype,
			const std::vector<Property> elementDOFs,
			const std::vector<Property> faceDOFs,
			const std::vector<Property> edgeDOFs,
			const std::vector<Property> pointDOFs,
			const std::vector<Property> midPointDOFs)
	: Physics(mesh, constraints, configuration, mtype, elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs) {};
	virtual ~TransientPhysics()
	{
		M.resize(_mesh.parts());
	}

	std::vector<SparseMatrix> M;
	std::vector<double> A;

protected:
	virtual void composeSubdomain(size_t subdomain) =0;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_TRANSIENT_ASSEMBLER_H_ */
