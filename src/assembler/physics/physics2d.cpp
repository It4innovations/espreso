
#include "physics2d.h"

#include "../../mesh/structures/mesh.h"
#include "../constraints/equalityconstraints.h"

using namespace espreso;

void Physics2D::prepareHybridTotalFETIWithCorners()
{
	prepare();
	_mesh->computePlaneCorners(1, true, true);
}

void Physics2D::prepareHybridTotalFETIWithKernels()
{
	prepare();
	_mesh->computeEdgesSharedByDomains();
}

