
#include "builder.utils.h"

namespace espreso {
namespace builder {

void fillNodes(ClusteredMesh &linked, Mesh &mesh)
{
//	size_t threads = info::env::OMP_NUM_THREADS;
//
//	info::mesh->nodes->size = _meshData.coordinates.size();
//	info::mesh->nodes->distribution = tarray<size_t>::distribute(threads, _meshData.coordinates.size());
//
//	info::mesh->nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(info::mesh->nodes->distribution, _meshData.nIDs));
//	info::mesh->nodes->coordinates = new serializededata<esint, Point >(1, tarray<Point>(info::mesh->nodes->distribution, _meshData.coordinates.cbegin(), _meshData.coordinates.cend()));
//
//	std::vector<size_t> rdistribution = info::mesh->nodes->distribution, rdatadistribution = info::mesh->nodes->distribution;
//	for (size_t t = 1; t < threads; t++) {
//		++rdistribution[t];
//		if (rdistribution[t] < _meshData._nrankdist.size()) {
//			rdatadistribution[t] = _meshData._nrankdist[rdistribution[t]];
//		} else {
//			rdatadistribution[t] = _meshData._nrankdist[rdistribution[threads] - 1];
//		}
//	}
//	++rdistribution[threads];
//	rdatadistribution[threads] = _meshData._nrankdist[rdistribution[threads] - 1];
//
//	info::mesh->nodes->ranks = new serializededata<esint, int>(tarray<esint>(rdistribution, _meshData._nrankdist), tarray<int>(rdatadistribution, _meshData._nranks));
//
//	info::mesh->boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES"));
//	info::mesh->boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, _meshData.nIDs.size()));
//	std::iota(info::mesh->boundaryRegions.back()->nodes->datatarray().begin(), info::mesh->boundaryRegions.back()->nodes->datatarray().end(), 0);
}

void fillElements(ClusteredMesh &linked, Mesh &mesh)
{

}

}
}
