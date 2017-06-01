
#include "physics3d.h"

#include "../../mesh/structures/mesh.h"
#include "../constraints/equalityconstraints.h"

using namespace espreso;

void Physics3D::prepareHybridTotalFETIWithCorners()
{
	prepareTotalFETI();
	_mesh->computeVolumeCorners(1, true, true, false);
}

void Physics3D::prepareHybridTotalFETIWithKernels()
{
	prepareTotalFETI();
	_mesh->computeFacesSharedByDomains();
}

void Physics3D::assembleB0FromCorners(const Step &step)
{
	EqualityConstraints::insertCornersGluingToB0(*_instance, _mesh->corners(), _nodesDOFsOffsets);
}

void Physics3D::assembleB0FromKernels(const Step &step)
{
	EqualityConstraints::insertKernelsGluingToB0(*_instance, _mesh->faces(), _mesh->nodes(), _nodesDOFsOffsets);
}



