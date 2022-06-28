
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
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

static serializededata<esint, esint>* outputOffset(const std::vector<DataDuplication> &duplication, const serializededata<esint, esint> *offset)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > odist(threads), odata(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		std::vector<esint> tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}
		auto dup = std::lower_bound(duplication.begin(), duplication.end(), DataDuplication{-1, *offset->datatarray().cbegin(t)});
		for (auto IDs = offset->datatarray().cbegin(t); IDs != offset->datatarray().cend(t); ++IDs) {
			while (dup != duplication.end() && dup->origin < *IDs) {
				++dup; // skip removed nodes
			}
			bool pushOrigin = true;
			while (dup != duplication.end() && dup->origin == *IDs) {
				if (dup->duplicate > dup->origin) {
					tdata.push_back(dup->origin);
					pushOrigin = false;
				}
				tdata.push_back(dup->duplicate);
				++dup;
			}
			if (pushOrigin) {
				tdata.push_back(*IDs);
			}
			tdist.push_back(tdata.size());
		}
		odist[t].swap(tdist);
		odata[t].swap(tdata);
	}
	utils::threadDistributionToFullDistribution(odist);
	return new serializededata<esint, esint>(odist, odata);
}

void fillNodes(const MergedNodes *nodes, OrderedRegions &regions, Mesh &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	mesh.nodes->size = nodes->offsets.size();
	mesh.nodes->distribution = tarray<size_t>::distribute(threads, nodes->offsets.size());

	mesh.nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(mesh.nodes->distribution, nodes->offsets.begin(), nodes->offsets.end()));
	mesh.nodes->coordinates = new serializededata<esint, Point >(1, tarray<Point>(mesh.nodes->distribution, nodes->coordinates.cbegin(), nodes->coordinates.cend()));
	mesh.nodes->outputOffset = outputOffset(nodes->duplication, mesh.nodes->IDs);

	mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES"));
	mesh.boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, mesh.nodes->size));
	std::iota(mesh.boundaryRegions.back()->nodes->datatarray().begin(), mesh.boundaryRegions.back()->nodes->datatarray().end(), 0);

	for (auto region = regions.nodes.begin(); region != regions.nodes.end(); ++region) {
		mesh.boundaryRegions.push_back(new BoundaryRegionStore(region->name));
		esint rsize = 0;
		for (size_t i = 0; i < nodes->offsets.size(); ++i) {
			if (region->start <= nodes->offsets[i] && nodes->offsets[i] < region->end) {
				++rsize;
			}
		}
		mesh.boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, rsize));
		for (size_t i = 0, j = 0; i < nodes->offsets.size(); ++i) {
			if (region->start <= nodes->offsets[i] && nodes->offsets[i] < region->end) {
				mesh.boundaryRegions.back()->nodes->datatarray()[j++] = i;
			}
		}
	}
}

void fillElements(const MergedElements *elements, OrderedRegions &regions, Mesh &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	size_t esize = 0;
	for (size_t e = 0; e < elements->etype.size(); ++e) {
		if (Mesh::element(elements->etype[e]).dimension == mesh.dimension) {
			++esize;
		}
	}

	esint eoffset = esize;
	esint totalSize = Communication::exscan(eoffset);
	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, esize);

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads), toffset(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	esize = 0;
	tedist.front().push_back(0);
	for (size_t e = 0, t = 0; e < elements->etype.size(); ++e) {
		if (Mesh::element(elements->etype[e]).dimension == mesh.dimension) {
			toffset[t].push_back(eoffset + esize); // = elements->offsets[e]
			for (esint n = 0; n < Mesh::element(elements->etype[e]).nodes; ++n) {
				tnodes[t].push_back(elements->enodes[elements->edist[e] + n]);
			}
			tedist[t].push_back(tnodes[t].size());
			epointers[t].push_back(&Mesh::edata[(int)elements->etype[e]]);

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
	mesh.elements->offset = new serializededata<esint, esint>(1, toffset);
	mesh.elements->nodes = new serializededata<esint, esint>(tedist, tnodes);
	mesh.elements->epointers = new serializededata<esint, Element*>(1, epointers);
	mesh.elements->outputOffset = outputOffset(elements->duplication, mesh.elements->offset);

	mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, esize));
	std::iota(mesh.elementsRegions.back()->elements->datatarray().begin(), mesh.elementsRegions.back()->elements->datatarray().end(), 0);

	for (auto region = regions.elements.begin(); region != regions.elements.end(); ++region) {
		esint rsize = 0, rnodes = 0;
		for (size_t i = 0; i < elements->offsets.size(); ++i) {
			if (region->start <= elements->offsets[i] && elements->offsets[i] < region->end) {
				++rsize;
				rnodes += Mesh::element(elements->etype[i]).nodes;
			}
		}
		if (region->dimension == mesh.dimension) {
			mesh.elementsRegions.push_back(new ElementsRegionStore(region->name));
			mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, rsize));
			for (size_t i = 0, j = 0; i < elements->offsets.size(); ++i) {
				if (region->start <= elements->offsets[i] && elements->offsets[i] < region->end) {
					mesh.elementsRegions.back()->elements->datatarray()[j++] = i;
				}
			}
		} else {
			std::vector<std::vector<esint> > tbedist(threads), tbnodes(threads);
			std::vector<std::vector<Element*> > tbepointers(threads);
			std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, rsize);
			std::vector<esint> eregiondist(rsize + 1);
			for (size_t i = 0, e = 0; i < elements->offsets.size(); ++i) {
				if (region->start <= elements->offsets[i] && elements->offsets[i] < region->end) {
					eregiondist[e + 1] = eregiondist[e] + Mesh::element(elements->etype[i]).nodes;
					++e;
				}
			}

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				tbedist[t].clear();
				if (t == 0) {
					tbedist[t].insert(tbedist[t].end(), eregiondist.begin() + edistribution[t], eregiondist.begin() + edistribution[t + 1] + 1);
				} else {
					tbedist[t].insert(tbedist[t].end(), eregiondist.begin() + edistribution[t] + 1, eregiondist.begin() + edistribution[t + 1] + 1);
				}

				tbnodes[t].resize(eregiondist[edistribution[t + 1]] - eregiondist[edistribution[t]]);
				tbepointers[t].resize(edistribution[t + 1] - edistribution[t]);
				for (size_t i = 0, e = 0, index = 0; i < elements->offsets.size(); ++i) {
					if (region->start <= elements->offsets[i] && elements->offsets[i] < region->end) {
						tbepointers[t][e] = &Mesh::edata[(int)elements->etype[i]];
						for (esint n = elements->edist[i]; n < elements->edist[i + 1]; ++n, ++index) {
							tbnodes[t][index] = elements->enodes[n];
						}
						++e;
					}
				}
			}

			mesh.boundaryRegions.push_back(new BoundaryRegionStore(region->name));
			mesh.boundaryRegions.back()->dimension = region->dimension;
			mesh.boundaryRegions.back()->originalDimension = region->dimension;
			mesh.boundaryRegions.back()->distribution.threads = edistribution;
			mesh.boundaryRegions.back()->elements = new serializededata<esint, esint>(tbedist, tbnodes);
			mesh.boundaryRegions.back()->epointers = new serializededata<esint, Element*>(1, tbepointers);
		}
	}
}

void fillSequentialMesh(const TemporalSequentialMesh<MergedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	fillNodes(prepared.nodes, regions, mesh);

	std::vector<size_t> rdistribution = mesh.nodes->distribution, rdatadistribution = mesh.nodes->distribution;
	++rdistribution[threads];

	mesh.nodes->ranks = new serializededata<esint, int>(tarray<esint>(rdistribution, 1UL, 0), tarray<int>(rdatadistribution, 1UL, 0));
	std::iota(mesh.nodes->ranks->boundarytarray().begin(), mesh.nodes->ranks->boundarytarray().end(), 0);

	fillElements(prepared.elements, regions, mesh);

	mesh.neighbors = { };
	mesh.neighborsWithMe = { 0 };

	utils::clearVectors(prepared.nodes->coordinates, prepared.nodes->offsets, prepared.nodes->duplication);
	utils::clearVectors(prepared.elements->offsets, prepared.elements->etype, prepared.elements->enodes, prepared.elements->edist, prepared.elements->duplication);
}

void fillMesh(const TemporalMesh<LinkedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh)
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

	fillNodes(prepared.nodes, regions, mesh);

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


	fillElements(prepared.elements, regions, mesh);

	mesh.neighbors = prepared.nodes->neighbors;
	mesh.neighborsWithMe.resize(mesh.neighbors.size() + 1);
	std::merge(mesh.neighbors.begin(), mesh.neighbors.end(), &info::mpi::rank, &info::mpi::rank + 1, mesh.neighborsWithMe.begin());

	utils::clearVectors(prepared.nodes->coordinates, prepared.nodes->offsets, prepared.nodes->duplication, prepared.nodes->rankData, prepared.nodes->rankDistribution, prepared.nodes->neighbors);
	utils::clearVectors(prepared.elements->offsets, prepared.elements->etype, prepared.elements->enodes, prepared.elements->edist, prepared.elements->duplication, prepared.elements->duplication);
}

}
}
