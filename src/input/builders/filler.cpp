
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "esinfo/envinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "wrappers/mpi/communication.h"

#include <numeric>

namespace espreso {
namespace builder {

void fillMesh(TemporalMesh<LinkedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	// nodes

	mesh.nodes->size = prepared.nodes->offsets.size();
	mesh.nodes->distribution = tarray<size_t>::distribute(threads, prepared.nodes->offsets.size());

	mesh.nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(mesh.nodes->distribution, prepared.nodes->offsets.begin(), prepared.nodes->offsets.end()));
	mesh.nodes->coordinates = new serializededata<esint, Point >(1, tarray<Point>(mesh.nodes->distribution, prepared.nodes->coordinates.cbegin(), prepared.nodes->coordinates.cend()));

	std::vector<size_t> rdistribution = mesh.nodes->distribution, rdatadistribution = mesh.nodes->distribution;
	for (size_t t = 1; t < threads; t++) {
		++rdistribution[t];
		if (rdistribution[t] < prepared.nodes->rankDistribution.size()) {
			rdatadistribution[t] = prepared.nodes->rankDistribution[rdistribution[t]];
		} else {
			rdatadistribution[t] = prepared.nodes->rankDistribution[rdistribution[threads] - 1];
		}
	}
	++rdistribution[threads];
	rdatadistribution[threads] = prepared.nodes->rankDistribution[rdistribution[threads] - 1];

	mesh.nodes->ranks = new serializededata<esint, int>(
			tarray<esint>(rdistribution, prepared.nodes->rankDistribution.begin(), prepared.nodes->rankDistribution.end()),
			tarray<int>(rdatadistribution, prepared.nodes->rankData.begin(), prepared.nodes->rankData.end()));

	mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES"));
	mesh.boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, mesh.nodes->size));
	std::iota(mesh.boundaryRegions.back()->nodes->datatarray().begin(), mesh.boundaryRegions.back()->nodes->datatarray().end(), 0);

	// elements

	size_t esize = 0;
	for (size_t e = 0; e < prepared.elements->etype.size(); ++e) {
		if (Mesh::element(prepared.elements->etype[e]).dimension == mesh.dimension) {
			++esize;
		}
	}

	esint eoffset = esize;
	esint totalSize = Communication::exscan(eoffset);
	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, esize);

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads), eIDs(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	esize = 0;
	tedist.front().push_back(0);
	for (size_t e = 0, t = 0; e < prepared.elements->etype.size(); ++e) {
		if (Mesh::element(prepared.elements->etype[e]).dimension == mesh.dimension) {
			eIDs[t].push_back(eoffset + esize); // = prepared.elements->offsets[e]
			for (esint n = 0; n < Mesh::element(prepared.elements->etype[e]).nodes; ++n) {
				tnodes[t].push_back(prepared.elements->enodes[prepared.elements->edist[e] + n]);
			}
			tedist[t].push_back(tnodes[t].size());
			epointers[t].push_back(&Mesh::edata[(int)prepared.elements->etype[e]]);

			++esize;
			if (esize == edistribution[t + 1]) {
				++t;
			}
		}
	}

	mesh.elements->distribution.process.offset = eoffset;
	mesh.elements->distribution.process.next = eoffset + esize;
	mesh.elements->distribution.process.size = esize;
	mesh.elements->distribution.process.totalSize = totalSize;
	mesh.elements->distribution.threads = edistribution;
	mesh.elements->IDs = new serializededata<esint, esint>(1, eIDs);
	mesh.elements->nodes = new serializededata<esint, esint>(tedist, tnodes);
	mesh.elements->epointers = new serializededata<esint, Element*>(1, epointers);

	mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, esize));
	std::iota(mesh.elementsRegions.back()->elements->datatarray().begin(), mesh.elementsRegions.back()->elements->datatarray().end(), 0);

	mesh.neighbors = prepared.nodes->neighbors;
	mesh.neighborsWithMe.resize(mesh.neighbors.size() + 1);
	std::merge(mesh.neighbors.begin(), mesh.neighbors.end(), &info::mpi::rank, &info::mpi::rank + 1, mesh.neighborsWithMe.begin());
}

}
}
