
#include "assembler.h"

#include "../../../mesh/structures/mesh.h"
#include "../../../basis/logging/logging.h"

#include "../../../solver/generic/SparseMatrix.h"

using namespace espreso;


void LinearPhysics::assembleStiffnessMatrices()
{
	ESINFO(PROGRESS2) << "Assemble matrices K, kernels, and RHS.";
	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		composeSubdomain(p);
		K[p].mtype = mtype;
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

