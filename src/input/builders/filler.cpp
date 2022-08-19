
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

void fillNodes(LinkedNodes &nodes, OrderedRegions &regions, Mesh &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	mesh.nodes->size = nodes.offsets.size();
	mesh.nodes->distribution = tarray<size_t>::distribute(threads, nodes.offsets.size());

	mesh.nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(mesh.nodes->distribution, nodes.offsets.begin(), nodes.offsets.end()));
	mesh.nodes->coordinates = new serializededata<esint, Point >(1, tarray<Point>(mesh.nodes->distribution, nodes.coordinates.cbegin(), nodes.coordinates.cend()));
	mesh.nodes->inputOffset = outputOffset(nodes.duplication, mesh.nodes->IDs);

	mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES"));
	mesh.boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, mesh.nodes->size));
	std::iota(mesh.boundaryRegions.back()->nodes->datatarray().begin(), mesh.boundaryRegions.back()->nodes->datatarray().end(), 0);

	for (auto region = regions.nodes.begin(); region != regions.nodes.end(); ++region) {
		mesh.boundaryRegions.push_back(new BoundaryRegionStore(region->name));
		esint rsize = 0;
		for (size_t i = 0; i < nodes.offsets.size(); ++i) {
			if (region->start <= nodes.offsets[i] && nodes.offsets[i] < region->end) {
				++rsize;
			}
		}
		mesh.boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, rsize));
		for (size_t i = 0, j = 0; i < nodes.offsets.size(); ++i) {
			if (region->start <= nodes.offsets[i] && nodes.offsets[i] < region->end) {
				mesh.boundaryRegions.back()->nodes->datatarray()[j++] = i;
			}
		}
	}

	std::vector<size_t> rdistribution = mesh.nodes->distribution, rdatadistribution = mesh.nodes->distribution;
	for (size_t t = 1; t < threads; t++) {
		++rdistribution[t];
		if (rdistribution[t] < nodes.rankDistribution.size()) {
			rdatadistribution[t] = nodes.rankDistribution[rdistribution[t]];
		} else {
			rdatadistribution[t] = nodes.rankDistribution[rdistribution[threads] - 1];
		}
	}
	++rdistribution[threads];
	rdatadistribution[threads] = nodes.rankDistribution[rdistribution[threads] - 1];

	mesh.nodes->ranks = new serializededata<esint, int>(
			tarray<esint>(rdistribution, nodes.rankDistribution.begin(), nodes.rankDistribution.end()),
			tarray<int>(rdatadistribution, nodes.rankData.begin(), nodes.rankData.end()));

	mesh.neighbors = nodes.neighbors;
	mesh.neighborsWithMe.resize(mesh.neighbors.size() + 1);
	std::merge(mesh.neighbors.begin(), mesh.neighbors.end(), &info::mpi::rank, &info::mpi::rank + 1, mesh.neighborsWithMe.begin());
}

void fillElements(MergedElements &elements, OrderedRegions &regions, Mesh &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	size_t esize = 0;
	for (size_t e = 0; e < elements.etype.size(); ++e) {
		if (Mesh::element(elements.etype[e]).dimension == mesh.dimension) {
			++esize;
		}
	}

	esint eoffset = esize;
	esint totalSize = Communication::exscan(eoffset);
	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, esize);

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads), toffset(threads), toutput(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	esize = 0;
	tedist.front().push_back(0);
	for (size_t e = 0, enodes = 0, t = 0; e < elements.etype.size(); enodes += Element::encode(elements.etype[e++]).nodes) {
		if (Mesh::element(elements.etype[e]).dimension == mesh.dimension) {
			toffset[t].push_back(eoffset + esize); // = elements->offsets[e]
			for (esint n = 0; n < Element::encode(elements.etype[e]).nodes; ++n) {
				tnodes[t].push_back(elements.enodes[enodes + n]);
			}
			tedist[t].push_back(tnodes[t].size());
			epointers[t].push_back(&Mesh::edata[(int)Element::encode(elements.etype[e]).code]);

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

	serializededata<esint, esint> out(1, tarray<esint>(edistribution, elements.offsets.begin(), elements.offsets.end())); // remove
	mesh.elements->inputOffset = outputOffset(elements.duplication, &out);

	mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, esize));
	std::iota(mesh.elementsRegions.back()->elements->datatarray().begin(), mesh.elementsRegions.back()->elements->datatarray().end(), 0);

	for (auto region = regions.elements.begin(); region != regions.elements.end(); ++region) {
		esint rsize = 0, rnodes = 0;
		for (size_t i = 0; i < elements.offsets.size(); ++i) {
			if (region->start <= elements.offsets[i] && elements.offsets[i] < region->end) {
				++rsize;
				rnodes += Element::encode(elements.etype[i]).nodes;
			}
		}
		if (region->dimension == mesh.dimension) {
			mesh.elementsRegions.push_back(new ElementsRegionStore(region->name));
			mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, rsize));
			for (size_t i = 0, j = 0; i < elements.offsets.size(); ++i) {
				if (region->start <= elements.offsets[i] && elements.offsets[i] < region->end) {
					mesh.elementsRegions.back()->elements->datatarray()[j++] = i;
				}
			}
		} else {
			// TODO: implement me
//			std::vector<esint> rnodes;
//			std::vector<Element*> repointers;
//			std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, rsize);
//			std::vector<esint> eregiondist(rsize + 1);
//
////			std::vector<std::vector<esint> > tbedist(threads), tbnodes(threads);
////			std::vector<std::vector<Element*> > tbepointers(threads);
////
////			std::vector<size_t> erdistribution = tarray<size_t>::distribute(threads, rsize);
////			std::vector<esint> eregiondist(rsize + 1);
////			std::vector<size_t> fullDistribution = { 0 };
////			for (size_t i = 0, e = 0; i < elements.offsets.size(); ++i) {
////				if (region->start <= elements.offsets[i] && elements.offsets[i] < region->end) {
////					eregiondist[e + 1] = eregiondist[e] + Mesh::element(elements.etype[i]).nodes;
////					if (e == 0) {
////						fullDistribution.front() = i;
////					}
////					++e;
////				}
////				if (erdistribution[fullDistribution.size()] == e) {
////					fullDistribution.push_back(i + 1);
////				}
////			}
////
////			#pragma omp parallel for
////			for (size_t t = 0; t < threads; t++) {
////				tbedist[t].clear();
////				if (t == 0) {
////					tbedist[t].insert(tbedist[t].end(), eregiondist.begin() + erdistribution[t], eregiondist.begin() + erdistribution[t + 1] + 1);
////				} else {
////					tbedist[t].insert(tbedist[t].end(), eregiondist.begin() + erdistribution[t] + 1, eregiondist.begin() + erdistribution[t + 1] + 1);
////				}
////
////				tbnodes[t].resize(eregiondist[erdistribution[t + 1]] - eregiondist[erdistribution[t]]);
////				tbepointers[t].resize(erdistribution[t + 1] - erdistribution[t]);
////				for (size_t i = fullDistribution[t], e = 0, index = 0; i < fullDistribution[t + 1]; ++i) {
////					if (region->start <= elements.offsets[i] && elements.offsets[i] < region->end) {
////						tbepointers[t][e++] = &Mesh::edata[(int)elements.etype[i]];
////						auto enodes = mesh.elements->nodes->cbegin() + i;
////						for (auto n = enodes->begin(); n != enodes->end(); ++n, ++index) {
////							tbnodes[t][index] = *n;
////						}
////					}
////				}
////			}
//
//			mesh.boundaryRegions.push_back(new BoundaryRegionStore(region->name));
//			mesh.boundaryRegions.back()->dimension = region->dimension;
//			mesh.boundaryRegions.back()->originalDimension = region->dimension;
//			mesh.boundaryRegions.back()->distribution.threads = rdistribution;
//			mesh.boundaryRegions.back()->elements = new serializededata<esint, esint>(tbedist, tbnodes);
//			mesh.boundaryRegions.back()->epointers = new serializededata<esint, Element*>(1, repointers);
//		}
		}
	}
}

}
}
