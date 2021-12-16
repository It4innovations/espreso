
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "esinfo/envinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "wrappers/mpi/communication.h"

#include <numeric>
//#include <fstream>

namespace espreso {
namespace builder {

void fillMesh(TemporalSequentialMesh<MergedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	// nodes

	mesh.nodes->size = prepared.nodes->offsets.size();
	mesh.nodes->distribution = tarray<size_t>::distribute(threads, prepared.nodes->offsets.size());

	mesh.nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(mesh.nodes->distribution, prepared.nodes->offsets.begin(), prepared.nodes->offsets.end()));
	mesh.nodes->coordinates = new serializededata<esint, Point >(1, tarray<Point>(mesh.nodes->distribution, prepared.nodes->coordinates.cbegin(), prepared.nodes->coordinates.cend()));

	std::vector<size_t> rdistribution = mesh.nodes->distribution, rdatadistribution = mesh.nodes->distribution;
	++rdistribution[threads];

	mesh.nodes->ranks = new serializededata<esint, int>(tarray<esint>(rdistribution, 1UL, 0), tarray<int>(rdatadistribution, 1UL, 0));
	std::iota(mesh.nodes->ranks->boundarytarray().begin(), mesh.nodes->ranks->boundarytarray().end(), 0);

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
	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, esize);

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads), eIDs(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	esize = 0;
	tedist.front().push_back(0);
	for (size_t e = 0, t = 0; e < prepared.elements->etype.size(); ++e) {
		if (Mesh::element(prepared.elements->etype[e]).dimension == mesh.dimension) {
			eIDs[t].push_back(esize); // = prepared.elements->offsets[e]
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

	mesh.elements->distribution.process.offset = 0;
	mesh.elements->distribution.process.next = esize;
	mesh.elements->distribution.process.size = esize;
	mesh.elements->distribution.process.totalSize = esize;
	mesh.elements->distribution.threads = edistribution;
	mesh.elements->IDs = new serializededata<esint, esint>(1, eIDs);
	mesh.elements->nodes = new serializededata<esint, esint>(tedist, tnodes);
	mesh.elements->epointers = new serializededata<esint, Element*>(1, epointers);

	mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, esize));
	std::iota(mesh.elementsRegions.back()->elements->datatarray().begin(), mesh.elementsRegions.back()->elements->datatarray().end(), 0);

	mesh.neighbors = { };
	mesh.neighborsWithMe = { 0 };

	prepared.clear();
}

void fillMesh(TemporalMesh<LinkedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh)
{
//	Communication::serialize([&] () {
//		for (size_t i = 0; i < prepared.nodes->neighbors.size(); ++i) {
//			std::ofstream os("ranks_" + std::to_string(std::min(info::mpi::rank, prepared.nodes->neighbors[i])) + "-" + std::to_string(std::max(info::mpi::rank, prepared.nodes->neighbors[i])) + "." + std::to_string(info::mpi::rank) + ".txt");
//			for (size_t n = 0; n < prepared.nodes->offsets.size(); ++n) {
//				for (esint r = prepared.nodes->rankDistribution[n]; r < prepared.nodes->rankDistribution[n + 1]; ++r) {
//					if (prepared.nodes->rankData[r] != info::mpi::rank && !std::binary_search(prepared.nodes->neighbors.begin(), prepared.nodes->neighbors.end(), prepared.nodes->rankData[r])) {
//						printf("invalid rank %d on process %d\n", prepared.nodes->rankData[r], info::mpi::rank);
//					}
//					if (prepared.nodes->rankData[r] == prepared.nodes->neighbors[i]) {
//						os << prepared.nodes->offsets[n] << "\n";
//					}
//				}
//			}
//		}
//
//		printf(" -- %d -- \n", info::mpi::rank);
//		for (size_t i = 0; i < prepared.nodes->offsets.size(); ++i) {
//			printf("%d: %.2f %.2f\n", prepared.nodes->offsets[i], prepared.nodes->coordinates[i].x, prepared.nodes->coordinates[i].y);
//		}
//
//		for (size_t i = 0; i < prepared.elements->offsets.size(); ++i) {
//			printf("%d:", prepared.elements->offsets[i]);
//			for (esint n = prepared.elements->edist[i]; n < prepared.elements->edist[i + 1]; ++n) {
//				printf(" %d", prepared.nodes->offsets[prepared.elements->enodes[n]]);
//			}
//			printf("\n");
//		}
//	});
//
//	if (!std::is_sorted(prepared.nodes->offsets.begin(), prepared.nodes->offsets.end())) {
//		printf("nodes on process %d are not sorted\n", info::mpi::rank);
//	}

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

	prepared.clear();
}

}
}
