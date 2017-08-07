
#include "physics3d.h"

#include "../../mesh/structures/mesh.h"
#include "../constraints/equalityconstraints.h"

using namespace espreso;

void Physics3D::prepareHybridTotalFETIWithCorners()
{
	prepare();
	_mesh->computeVolumeCorners(1, true, true, false);
}

void Physics3D::prepareHybridTotalFETIWithKernels()
{
	prepare();
	_mesh->computeFacesSharedByDomains();
}



