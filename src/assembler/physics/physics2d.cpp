
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

void Physics2D::assembleB0FromCorners()
{
	EqualityConstraints::insertCornersGluingToB0(*_instance, _mesh->corners(), _nodesDOFsOffsets);
}

void Physics2D::assembleB0FromKernels(const std::vector<SparseMatrix> &kernels)
{
	EqualityConstraints::insertKernelsGluingToB0(*_instance, _mesh->edges(), _mesh->nodes(), _nodesDOFsOffsets, kernels);
}



