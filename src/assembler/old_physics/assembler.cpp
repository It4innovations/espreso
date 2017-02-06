
#include "../old_physics/assembler.h"

#include "../../mesh/structures/mesh.h"
#include "../../basis/utilities/utils.h"

#include "../../solver/generic/SparseMatrix.h"

#include "../../config/solverespreso.h"

using namespace espreso;

OldPhysics::OldPhysics(Mesh &mesh,
		Constraints &constraints,
		const ESPRESOSolver &configuration,
		MatrixType mtype,
		const std::vector<Property> elementDOFs,
		const std::vector<Property> faceDOFs,
		const std::vector<Property> edgeDOFs,
		const std::vector<Property> pointDOFs,
		const std::vector<Property> midPointDOFs):
elementDOFs(elementDOFs),
faceDOFs(faceDOFs),
edgeDOFs(edgeDOFs),
pointDOFs(pointDOFs),
midPointDOFs(midPointDOFs),
matrixSize(mesh.parts()),
mtype(mtype),
_mesh(mesh),
_constraints(constraints),
_solverConfiguration(configuration)
{
	K.resize(_mesh.parts());
	R1.resize(_mesh.parts());
	R2.resize(_mesh.parts());
	RegMat.resize(_mesh.parts());
	D.resize(_mesh.parts());
	f.resize(_mesh.parts());
	singularK.resize(_mesh.parts(), true);
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		ESTEST(MANDATORY) << "Do not use HYBRID FETI for clusters with 1 domain." << (_mesh.parts() > 1 ? TEST_PASSED : TEST_FAILED);
	}
}

OldPhysics::~OldPhysics()
{

}

void OldPhysics::saveStiffnessMatrices()
{
	ESINFO(PROGRESS2) << "Save matrices K and RHS.";
	for (size_t p = 0; p < K.size(); p++) {
		std::ofstream osK(Logging::prepareFile(p, "K").c_str());
		osK << K[p];
		osK.close();
	}

	for (size_t p = 0; p < D.size(); p++) {
		std::ofstream osD(Logging::prepareFile(p, "D").c_str());
		osD << D[p];
		osD.close();
	}

	for (size_t p = 0; p < f.size(); p++) {
		std::ofstream osF(Logging::prepareFile(p, "f").c_str());
		osF << f[p];
		osF.close();
	}
}

void OldPhysics::saveKernelMatrices()
{
	for (size_t p = 0; p < R1.size(); p++) {
		std::ofstream osR(Logging::prepareFile(p, "R1").c_str());
		SparseMatrix tmpR = R1[p];
		tmpR.ConvertDenseToCSR(0);
		osR << tmpR;
		osR.close();
	}

	for (size_t p = 0; p < R2.size(); p++) {
		std::ofstream osR(Logging::prepareFile(p, "R2").c_str());
		SparseMatrix tmpR = R2[p];
		tmpR.ConvertDenseToCSR(0);
		osR << tmpR;
		osR.close();
	}

	for (size_t p = 0; p < RegMat.size(); p++) {
		std::ofstream osR(Logging::prepareFile(p, "RegMat").c_str());
		osR << RegMat[p];
		osR.close();
	}
}


