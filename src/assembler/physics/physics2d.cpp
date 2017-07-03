
#include "physics2d.h"

#include "../../mesh/structures/mesh.h"
#include "../constraints/equalityconstraints.h"

using namespace espreso;

void Physics2D::prepareHybridTotalFETIWithCorners()
{
	prepareTotalFETI();
	_mesh->computePlaneCorners(1, true, true);
}

void Physics2D::prepareHybridTotalFETIWithKernels()
{
	prepareTotalFETI();
	_mesh->computeEdgesSharedByDomains();
}

void Physics2D::assembleB0FromCorners(const Step &step)
{
	EqualityConstraints::insertCornersGluingToB0(*_instance, _mesh->corners(), _nodesDOFsOffsets);
}

void Physics2D::assembleB0FromKernels(const Step &step, const std::vector<SparseMatrix> &kernels)
{
	EqualityConstraints::insertKernelsGluingToB0(*_instance, _mesh->edges(), _mesh->nodes(), _nodesDOFsOffsets, kernels);
}



