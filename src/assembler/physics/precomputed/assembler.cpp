
#include "assembler.h"

#include "../../../basis/logging/logging.h"

#include "../../../solver/generic/SparseMatrix.h"
#include "../../../mesh/structures/mesh.h"

using namespace espreso;

PrecomputedPhysics::PrecomputedPhysics(
		APIMesh &mesh,
		Constraints &constraints,
		const ESPRESOSolver &configuration,
		MatrixType mtype,
		const std::vector<Property> elementDOFs,
		const std::vector<Property> faceDOFs,
		const std::vector<Property> edgeDOFs,
		const std::vector<Property> pointDOFs,
		const std::vector<Property> midPointDOFs,
		double *rhs, size_t rhs_size)
: Physics(mesh, constraints, configuration, mtype, elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs), _apimesh(mesh), _rhs(rhs), _rhs_size(rhs_size) {};

void PrecomputedPhysics::assembleStiffnessMatrices()
{
	ESINFO(PROGRESS2) << "Assemble matrices K and RHS.";
	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		composeSubdomain(p);
		K[p].mtype = mtype;
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

void PrecomputedPhysics::saveMeshProperties(store::Store &store)
{
	ESINFO(GLOBAL_ERROR) << "It is not possible to save mesh through API";
}

void PrecomputedPhysics::saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results)
{
	ESINFO(GLOBAL_ERROR) << "It is not possible to save results through API";
}




