
#include "meshpreprocessing.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "output/visualization/debug.h"

#include "wrappers/mpi/communication.h"
#include "wrappers/kahip/w.kahip.h"
#include "wrappers/metis/w.metis.h"
#include "wrappers/parmetis/w.parmetis.h"
#include "wrappers/ptscotch/w.ptscotch.h"
#include "wrappers/scotch/w.scotch.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <queue>

namespace espreso {
namespace mesh {

void computeNodesDuplication(NodeStore *nodes, std::vector<int> &neighborsWithMe)
{
	profiler::syncstart("compute_nodes_duplication");

	if (nodes->ranks) {
		delete nodes->ranks;
	}

	std::vector<esint> nids(nodes->IDs->datatarray().begin(), nodes->IDs->datatarray().end());

	std::vector<std::vector<esint> > sBuffer(neighborsWithMe.size(), nids), rBuffer(neighborsWithMe.size());
	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighborsWithMe)) {
		eslog::internalFailure("cannot exchange nodes ids.\n");
	}

	int threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > rdist(threads);
	std::vector<std::vector<int> >rdata(threads);
	rdist.front().push_back(0);
	std::vector<size_t> offset(neighborsWithMe.size());
	for (size_t n = 0; n < nids.size(); ++n) {
		for (size_t r = 0; r < neighborsWithMe.size(); ++r) {
			while (offset[r] < rBuffer[r].size() && rBuffer[r][offset[r]] < nids[n]) { ++offset[r]; }
			if (offset[r] < rBuffer[r].size() && rBuffer[r][offset[r]] == nids[n]) {
				rdata.front().push_back(neighborsWithMe[r]);
			}
		}
		rdist.front().push_back(rdata.front().size());
	}

	serializededata<esint, int>::balance(rdist, rdata);
	nodes->ranks = new serializededata<esint, int>(rdist, rdata);

	profiler::syncend("compute_nodes_duplication");
	eslog::checkpointln("MESH: NODES DUPLICATION COMPUTED");
}

void linkNodesAndElements(ElementStore *elements, NodeStore *nodes, const std::vector<int> &neighbors)
{
	profiler::syncstart("link_nodes_and_elements");

	std::vector<esint> n2eDist(1 + nodes->size + 1), _n2eData;

	auto nit = elements->nodes->cbegin();
	for (size_t e = 0; e < elements->epointers->datatarray().size(); ++e, ++nit) {
		PolyElement poly(Element::decode(elements->epointers->datatarray()[e]->code, nit->size()), nit->begin());
		for (auto n = nit->begin(); n != nit->end(); ++n) {
			if (poly.isNode(n - nit->begin())) {
				++n2eDist[*n + 1];
			}
		}
	}

	std::vector<esint> _n2eDist(n2eDist.begin(), n2eDist.end());
	utils::sizesToOffsets(_n2eDist);
	_n2eData.resize(_n2eDist.back(), -1);

	nit = elements->nodes->cbegin();
	for (size_t e = 0; e < elements->epointers->datatarray().size(); ++e, ++nit) {
		PolyElement poly(Element::decode(elements->epointers->datatarray()[e]->code, nit->size()), nit->begin());
		for (auto n = nit->begin(); n != nit->end(); ++n) {
			if (poly.isNode(n - nit->begin())) {
				_n2eData[_n2eDist[*n + 1]++] = elements->offset->datatarray()[e];
			}
		}
	}

	std::vector<std::vector<esint> > sBuffer(neighbors.size()), rBuffer(neighbors.size());
	auto ranks = nodes->ranks->cbegin();
	for (esint n = 0; n < nodes->size; ++n, ++ranks) {
		int r = 0;
		for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
			if (*rank != info::mpi::rank) {
				while (neighbors[r] < *rank) ++r;
				sBuffer[r].push_back(-nodes->IDs->datatarray()[n]); // nIDs are negative
				for (esint i = _n2eDist[n]; i < _n2eDist[n + 1]; ++i) {
					// polygons can put more elements to the same node
					if (i == _n2eDist[n] || _n2eData[i - 1] != _n2eData[i]) {
						sBuffer[r].push_back(_n2eData[i]);
					}
				}
			}
		}
		for (esint i = _n2eDist[n] + 1; i < _n2eDist[n + 1]; ++i) {
			if (_n2eData[i - 1] == _n2eData[i]) {
				--n2eDist[n + 1];
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighbors)) {
		eslog::internalFailure("addLinkFromTo - exchangeUnknownSize.\n");
	}
	utils::clearVector(sBuffer);

	for (size_t r = 0; r < rBuffer.size(); ++r) {
		auto it = rBuffer[r].begin();
		while (it != rBuffer[r].end()) {
			esint n = std::lower_bound(nodes->IDs->datatarray().begin(), nodes->IDs->datatarray().end(), -*it) - nodes->IDs->datatarray().begin();
			*it++ = -n;
			while (it != rBuffer[r].end() && 0 <= *it) {
				++n2eDist[n + 1];
				++it;
			}
		}
	}

	utils::sizesToOffsets(n2eDist);
	ivector<esint> n2eData(n2eDist.back());
	size_t neigh = 0;

	auto insertNeigh = [&] () {
		auto it = rBuffer[neigh].begin();
		while (it != rBuffer[neigh].end()) {
			esint n = -*it++;
			while (it != rBuffer[neigh].end() && 0 <= *it) {
				n2eData[n2eDist[n + 1]++] = *it;
				++it;
			}
		}
	};
	while (neigh < neighbors.size() && neighbors[neigh] < info::mpi::rank) {
		insertNeigh();
		++neigh;
	}
	for (esint n = 0; n < nodes->size; ++n) {
		for (esint i = _n2eDist[n]; i < _n2eDist[n + 1]; ++i) {
			if (i == _n2eDist[n] || _n2eData[i - 1] != _n2eData[i]) {
				n2eData[n2eDist[n + 1]++] = _n2eData[i];
			}
		}
	}
	while (neigh < neighbors.size()) {
		insertNeigh();
		++neigh;
	}
	n2eDist.pop_back();

	std::vector<size_t> bdist = nodes->distribution, ddist = nodes->distribution;
	for (size_t i = 1; i < nodes->distribution.size(); ++i) {
		ddist[i] = n2eDist[bdist[i]];
		bdist[i] += 1;
	}
	nodes->elements = new serializededata<esint, esint>(tarray<esint>(bdist, n2eDist), tarray<esint>(ddist, n2eData.begin(), n2eData.end()));

	profiler::synccheckpoint("rbuffer");
	profiler::syncend("link_nodes_and_elements");
	eslog::checkpointln("MESH: NODES AND ELEMENTS LINKED");
}

static void _computeElementsNeighbors(NodeStore *nodes, ElementStore *elements, bool faces)
{
	profiler::syncstart("compute_elements_neighbors");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > dualDistribution(threads);
	std::vector<std::vector<esint> > dualData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto enodes = elements->nodes->cbegin(t);
		const auto &edist = elements->offset->datatarray().distribution();

		std::vector<esint> tdist, tdata, intersection;
		tdist.reserve(edist[t + 1] - edist[t] + 1);
		tdata.reserve(tdist.capacity() * (faces ? 6 : 48));

		if (t == 0) {
			tdist.push_back(0);
		}

		auto init = [&intersection, &nodes] (esint element, esint node) {
			auto eit = nodes->elements->cbegin() + node;
			intersection.clear();
			for (auto n = eit->begin(); n != eit->end(); ++n) {
				if (*n != element) {
					intersection.push_back(*n);
				}
			}
		};

		auto merge = [&intersection, &nodes] (esint node) {
			auto eit = nodes->elements->cbegin() + node;
			auto it1 = intersection.begin();
			auto it2 = eit->begin();
			auto last = intersection.begin();
			while (it1 != intersection.end()) {
				while (it2 != eit->end() && *it2 < *it1) {
					++it2;
				}
				if (it2 == eit->end()) {
					break;
				}
				if (*it1 == *it2) {
					*last++ = *it1++;
				} else {
					it1++;
				}
			}
			intersection.resize(last - intersection.begin());
		};

		auto insert = [&] () {
			if (!faces) {
				tdata.push_back(intersection.size());
				tdata.insert(tdata.end(), intersection.begin(), intersection.end());
			} else {
				tdata.push_back(intersection.size() ? intersection.front() : -1);
				if (intersection.size() > 1) {
					eslog::error("Input error: a face shared by 3 elements found.\n");
				}
			}
		};

		for (size_t e = edist[t]; e < edist[t + 1]; ++e, ++enodes) {
			switch (Element::encode(elements->epointers->datatarray()[e]->code).code) {
			case Element::CODE::POLYGON:
				if (!faces) {
					eslog::internalFailure("implement polygon edge neighbors.\n");
				}
				break;
			case Element::CODE::POLYHEDRON:
				if (faces) {
					for (esint face = 0, i = 1; face < enodes->front(); ++face, i += enodes->at(i) + 1) {
						init(elements->offset->datatarray()[e], enodes->at(i + 1));
						for (esint n = 1; n < enodes->at(i); ++n) {
							merge(enodes->at(i + 1 + n));
						}
						insert();
					}
				} else {
					eslog::internalFailure("implement polyhedron edge neighbors.\n");
				}
				break;
			default: {
				serializededata<int, int> *interfaces = faces ? elements->epointers->datatarray()[e]->faces : elements->epointers->datatarray()[e]->edges;
				for (auto interface = interfaces->begin(); interface != interfaces->end(); ++interface) {
					init(elements->offset->datatarray()[e], enodes->at(*interface->begin()));
					for (auto n = interface->begin() + 1; n != interface->end() && intersection.size(); ++n) {
						merge(enodes->at(*n));
					}
					insert();
				}
			}
			}
			tdist.push_back(tdata.size());
		}

		dualDistribution[t].swap(tdist);
		dualData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dualDistribution);

	if (faces) {
		elements->faceNeighbors = new serializededata<esint, esint>(dualDistribution, dualData);
	} else {
		elements->edgeNeighbors = new serializededata<esint, esint>(dualDistribution, dualData);
	}

	profiler::syncend("compute_elements_neighbors");
	eslog::checkpointln("MESH: ELEMENTS NEIGHBOURS COMPUTED");
}

void computeElementsFaceNeighbors(NodeStore *nodes, ElementStore *elements)
{
	_computeElementsNeighbors(nodes, elements, true);
	DebugOutput::faceNeighbors();
}

void computeElementsEdgeNeighbors(NodeStore *nodes, ElementStore *elements)
{
	_computeElementsNeighbors(nodes, elements, false);
}

void computeElementsCenters(const NodeStore *nodes, ElementStore *elements)
{
	if (elements->centers) {
		return;
	}
	profiler::syncstart("compute_element_centers");
	int threads = info::env::OMP_NUM_THREADS;

	elements->centers = new serializededata<esint, Point>(1, elements->epointers->datatarray().distribution());

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		auto center = elements->centers->datatarray().begin(t);
		for (auto e = elements->nodes->cbegin(t); e != elements->nodes->cend(t); ++e, ++center) {
			for (auto n = e->begin(); n != e->end(); ++n) {
				*center += nodes->coordinates->datatarray()[*n];
			}
			*center /= e->size();
		}
	}
	profiler::syncend("compute_element_centers");
	eslog::checkpointln("MESH: ELEMENTS CENTERS COMPUTED");
}

esint getSFCDecomposition(const ElementStore *elements, const NodeStore *nodes, std::vector<esint> &partition)
{
	profiler::syncstart("get_sfc_decomposition");
	int threads = info::env::OMP_NUM_THREADS;

	HilbertCurve<double> sfc(info::mesh->dimension, SFCDEPTH, nodes->coordinates->datatarray().size(), nodes->coordinates->datatarray().data());

	std::vector<esint, initless_allocator<esint> > buckets(elements->epointers->datatarray().size());
	std::vector<esint, initless_allocator<esint> > permutation(elements->epointers->datatarray().size());
	std::vector<esint, initless_allocator<esint> > borders;

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		auto center = elements->centers->datatarray().cbegin(t);
		for (size_t e = elements->epointers->datatarray().distribution()[t]; e != elements->epointers->datatarray().distribution()[t + 1]; ++e, ++center) {
			buckets[e] = sfc.getBucket(*center);
		}
	}

	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		return buckets[i] < buckets[j];
	});
	profiler::synccheckpoint("sfc_buckets");

	if (!Communication::computeSplitters(buckets, permutation, borders, elements->distribution.process.totalSize, sfc.buckets(sfc.depth))) {
		eslog::internalFailure("cannot compute splitters.\n");
	}
	profiler::synccheckpoint("compute_splitters");

	if (elements->epointers->datatarray().size()) {
		auto border = std::lower_bound(borders.begin(), borders.end(), buckets[permutation[0]] + 1);
		for (size_t e = 0; e != elements->epointers->datatarray().size(); ++e) {
			while (*border < buckets[permutation[e]]) {
				++border;
			}
			partition[permutation[e]] = border - borders.begin() - 1;
		}
	}

	profiler::synccheckpoint("fill_partition");
	profiler::syncend("get_sfc_decomposition");
	return 0;

	std::vector<Point> dcenters(info::mpi::size), sumcenters(info::mpi::size);
	std::vector<esint> dsize(info::mpi::size), sumsize(info::mpi::size);

	for (size_t e = 0; e < elements->epointers->datatarray().size(); ++e) {
		dcenters[partition[e]] += elements->centers->datatarray()[e];
		++dsize[partition[e]];
	}

	Communication::allReduce(dcenters.data(), sumcenters.data(), 3 * info::mpi::size, MPI_DOUBLE, MPI_SUM);
	Communication::allReduce(dsize.data(), sumsize.data(), info::mpi::size, MPITools::getType<esint>().mpitype, MPI_SUM);

	for (int r = 0; r < info::mpi::size; ++r) {
		dcenters[r] = sumcenters[r] / sumsize[r];
	}

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		for (size_t e = elements->epointers->datatarray().distribution()[t]; e < elements->epointers->datatarray().distribution()[t + 1]; ++e) {
			Point &center = elements->centers->datatarray()[e];
			for (int r = 0; r < info::mpi::size; ++r) {
				if ((dcenters[r] - center).length() < (dcenters[partition[e]] - center).length()) {
					partition[e] = r;
				}
			}
		}
	}

	return 0; // edge cut is not computed
}

esint getVolumeDecomposition(const ElementStore *elements, const NodeStore *nodes, std::vector<esint> &partition)
{
	profiler::syncstart("get_sfc_decomposition");
	int threads = info::env::OMP_NUM_THREADS;

	HilbertCurve<double> sfc(info::mesh->dimension, SFCDEPTH, nodes->coordinates->datatarray().size(), nodes->coordinates->datatarray().data());

	std::vector<esint, initless_allocator<esint> > buckets(elements->epointers->datatarray().size());

	float vsum = 0;
	std::vector<float> volume(elements->epointers->datatarray().size());

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		size_t e = elements->epointers->datatarray().distribution()[t];
		for (auto enodes = elements->nodes->begin(t); enodes != elements->nodes->end(t); ++e, ++enodes) {
			buckets[e] = sfc.getBucket(nodes->coordinates->datatarray()[enodes->back()]);
			Point min = nodes->coordinates->datatarray()[enodes->back()], max = min;
			PolyElement poly(Element::decode(elements->epointers->datatarray()[e]->code, enodes->size()), enodes->data());
			for (size_t n = 0; n < enodes->size(); ++n) {
				if (poly.isNode(n)) {
					nodes->coordinates->datatarray()[enodes->at(n)].minmax(min, max);
				}
			}
			auto box = max - min;
			volume[e] = box.x * box.y * box.z;
			vsum += volume[e];
		}
	}

	struct __transfer__ {
		int from, to;
		float volume;
		int elements = 0;
	};

	struct __neigh__ {
		int rank;
		float volume;
	};

	std::vector<float> vdist = Communication::getDistribution(vsum);
	float vavg = vdist.back() / info::mpi::size;

	std::vector<std::vector<int> > nrank(info::mpi::size);
	std::vector<__transfer__> tr;
	std::queue<__neigh__> under, over;
	for (int r = 0; r < info::mpi::size; ++r) {
		if (vdist[r + 1] - vdist[r] < vavg) {
			under.push(__neigh__{ r, vavg - (vdist[r + 1] - vdist[r]) });
		} else {
			over.push(__neigh__{ r, (vdist[r + 1] - vdist[r]) - vavg });
		}
		while (!under.empty() && !over.empty()) {
			if (over.front().rank == info::mpi::rank || under.front().rank == info::mpi::rank) {
				tr.push_back(__transfer__{ over.front().rank, under.front().rank, std::min(under.front().volume, over.front().volume) });
			}
			if (under.front().volume < over.front().volume) {
				over.front().volume -= under.front().volume;
				under.pop();
			} else {
				under.front().volume -= over.front().volume;
				over.pop();
			}
		}
	}

	std::vector<esint, initless_allocator<esint> > permutation(elements->epointers->datatarray().size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		return volume[i] > volume[j];
	});

	float vprefix = 0;
	size_t prefix = 0;
	while (prefix < volume.size() && vprefix < vsum - vavg) {
		vprefix += volume[permutation[prefix++]];
	}
	if (prefix == 0) { // receivers exchange last SFC elements
		prefix = permutation.size();
	}
	std::sort(permutation.begin(), permutation.begin() + prefix, [&] (esint i, esint j) {
		return buckets[i] < buckets[j];
	});

	prefix = 0;
	std::vector<std::vector<int> > sBuffer(tr.size()), rBuffer(tr.size());
	std::vector<int> neighbors;
	for (size_t i = 0; i < tr.size(); ++i) {
		if (tr[i].from == info::mpi::rank) {
			vprefix = 0;
			while (prefix < volume.size() && vprefix < tr[i].volume) {
				vprefix += volume[permutation[prefix++]];
				++tr[i].elements;
			}
			sBuffer[neighbors.size()].push_back(tr[i].elements);
			neighbors.push_back(tr[i].to);
		} else {
			rBuffer[neighbors.size()].push_back(0);
			neighbors.push_back(tr[i].from);
		}
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, neighbors)) {
		eslog::error("cannot exchange volume-elements info.\n");
	}

	for (size_t i = 0; i < tr.size(); ++i) {
		if (tr[i].to == info::mpi::rank) {
			tr[i].elements = rBuffer[i].front();
		}
	}

	partition.clear();
	partition.resize(elements->epointers->datatarray().size(), info::mpi::rank);
	prefix = 0;
	for (size_t i = 0; i < tr.size(); ++i) {
		if (tr[i].from == info::mpi::rank) {
			for (int e = 0; e < tr[i].elements; ++e, ++prefix) {
				partition[permutation[prefix]] = tr[i].to;
			}
		} else {
			for (size_t e = permutation.size() - tr[i].elements; e < permutation.size(); ++e) {
				partition[permutation[e]] = tr[i].from;
			}
		}
	}

	profiler::syncend("get_sfc_decomposition");
	return 0;
}

void extractMeshDual(const ElementStore *elements, const NodeStore *nodes, std::vector<esint> &eframes, std::vector<esint> &eneighbors)
{
	// Disable due to horrible scalability
//	if (elements->centers == NULL) {
//		computeElementsCenters();
//	}

//	bool separateRegions = info::ecf->input.decomposition.separate_regions;
//	bool separateMaterials = info::ecf->input.decomposition.separate_materials;
//	bool separateEtypes = info::ecf->input.decomposition.separate_etypes;
//	esint eoffset = elements->distribution.process.offset;

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<esint> dDistribution(elements->epointers->datatarray().size() + 1);
	std::vector<std::vector<esint> > dData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;
//		int mat1 = 0, mat2 = 0, reg = 0, etype1 = 0, etype2 = 0;
//		int rsize = bitMastSize(elements->regions->datatarray().size() / elements->distribution.process.size);
//		esint hindex = 0;

		auto neighs = elements->faceNeighbors->cbegin(t);
		for (size_t e = elements->epointers->datatarray().distribution()[t]; e < elements->epointers->datatarray().distribution()[t + 1]; ++e, ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1) {
					tdata.push_back(*n);
//					if ((separateRegions || separateMaterials || separateEtypes) && (*n < eoffset || eoffset + elements->distribution.process.size <= *n)) {
//						hindex = std::lower_bound(halo->IDs->datatarray().begin(), halo->IDs->datatarray().end(), *n) - halo->IDs->datatarray().begin();
//					}
//					if (separateMaterials) {
//						mat1 = elements->material->datatarray()[e];
//						if (*n < eoffset || eoffset + elements->distribution.process.size <= *n) {
//							mat2 = halo->material->datatarray()[hindex];
//						} else {
//							mat2 = elements->material->datatarray()[*n - eoffset];
//						}
//					}
//					if (separateRegions) {
//						if (*n < eoffset || eoffset + elements->distribution.process.size <= *n) {
//							reg = memcmp(elements->regions->datatarray().data() + e * rsize, halo->regions->datatarray().data() + hindex * rsize, sizeof(esint) * rsize);
//						} else {
//							reg = memcmp(elements->regions->datatarray().data() + e * rsize, elements->regions->datatarray().data() + (*n - eoffset) * rsize, sizeof(esint) * rsize);
//						}
//					}
//					if (separateEtypes) {
//						etype1 = (int)elements->epointers->datatarray()[e]->type;
//						if (*n < eoffset || eoffset + elements->distribution.process.size <= *n) {
//							etype2 = (int)halo->epointers->datatarray()[hindex]->type;
//						} else {
//							etype2 = (int)elements->epointers->datatarray()[*n - eoffset]->type;
//						}
//					}
//
//					if (mat1 == mat2 && !reg && etype1 == etype2) {
//						tdata.push_back(*n);
//					}
				}
			}
			dDistribution[e + 1] = tdata.size();
		}

		dData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dDistribution, elements->epointers->datatarray().distribution());
	for (size_t t = 1; t < threads; t++) {
		dData[0].insert(dData[0].end(), dData[t].begin(), dData[t].end());
	}

	eframes.swap(dDistribution);
	eneighbors.swap(dData[0]);

	eslog::checkpointln("MESH: DUAL GRAPH EXTRACTED");
}

esint callParallelDecomposer(const ElementStore *elements, const NodeStore *nodes, std::vector<esint> &partition)
{
	esint edgecut = 0;

	if (
			(info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::NONE) ||
			(info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::METIS && !METIS::islinked()) ||
			(info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::PARMETIS && !ParMETIS::islinked()) ||
			(info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::PTSCOTCH && !PTScotch::islinked())
			) {

		eslog::checkpointln("MESH: RECLUSTERIZED SKIPPED");
		return edgecut;
	}

	if (info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE) {
		edgecut = getSFCDecomposition(elements, nodes, partition);
		eslog::checkpointln("MESH: RECLUSTERIZED ACCORDING SFC");
		return edgecut;
	}

	if (info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::VOLUME) {
		edgecut = getVolumeDecomposition(elements, nodes, partition);
		eslog::checkpointln("MESH: RECLUSTERIZED ACCORDING VOLUME");
		return edgecut;
	}

	std::vector<esint> eframes, eneighbors;
	extractMeshDual(elements, nodes, eframes, eneighbors);

	DebugOutput::meshDual(eframes, eneighbors);

	profiler::syncstart("call_parallel_decomposer");
	if (info::mpi::size <= info::ecf->input.third_party_scalability_limit && info::ecf->input.decomposition.parallel_decomposer != DecompositionConfiguration::ParallelDecomposer::METIS) {
		std::vector<esint> edistribution = Communication::getDistribution<esint>(partition.size(), &MPITools::subset->within);

		switch (info::ecf->input.decomposition.parallel_decomposer) {
		case DecompositionConfiguration::ParallelDecomposer::NONE: break;
		case DecompositionConfiguration::ParallelDecomposer::METIS: break; // never accessed
		case DecompositionConfiguration::ParallelDecomposer::PARMETIS:
			edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
					MPITools::subset->within, edistribution.data(), eframes.data(), eneighbors.data(),
					0, NULL, 0, NULL, NULL, partition.data());

			profiler::synccheckpoint("parmetis");
			eslog::checkpointln("MESH: RECLUSTERIZED BY PARMETIS");

			if (info::ecf->input.decomposition.parmetis_options.refinement) {
				esint prev = 2 * edgecut;
				while (1.01 * edgecut < prev) {
					prev = edgecut;
					edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_RefineKway,
							MPITools::subset->within, edistribution.data(), eframes.data(), eneighbors.data(),
							0, NULL, 0, NULL, NULL, partition.data());
				}
				profiler::synccheckpoint("parmetis");
				eslog::checkpointln("MESH: CLUSTERIZATION REFINED BY PARMETIS");
			}
			break;
		case DecompositionConfiguration::ParallelDecomposer::PTSCOTCH:
			edgecut = PTScotch::call(
					MPITools::subset->within, edistribution.data(), eframes.data(), eneighbors.data(),
					0, NULL, 0, NULL, NULL, partition.data());

			profiler::synccheckpoint("ptscotch");
			eslog::checkpointln("MESH: RECLUSTERIZED BY PTSCOTCH");
			break;
		case DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE: break; // never accessed
		case DecompositionConfiguration::ParallelDecomposer::VOLUME: break; // never accessed
		}
	} else {
		MPIType type = MPITools::getType<esint>();
		std::vector<esint> gframes, gneighbors, gpartition, edistribution;
		std::vector<size_t> offsets;

		if (info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::METIS) {
			Communication::gatherUnknownSize(eframes, gframes, &MPITools::singleton->across);
			Communication::gatherUnknownSize(eneighbors, gneighbors, &MPITools::singleton->across);
			Communication::gatherUnknownSize(partition, gpartition, offsets, &MPITools::singleton->across);
		} else {
			Communication::gatherUnknownSize(eframes, gframes, &MPITools::subset->across);
			Communication::gatherUnknownSize(eneighbors, gneighbors, &MPITools::subset->across);
			Communication::gatherUnknownSize(partition, gpartition, offsets, &MPITools::subset->across);
		}

		profiler::synccheckpoint("gather_dual");
		eslog::checkpointln("MESH: MPI PROCESSES REDUCED");

		auto fixframes = [&] () {
			size_t size = gframes.size();
			for (size_t i = 1, j = 1, offset = 0; j < gframes.size(); i++, j++) {
				while (gframes[j] == 0) {
					offset += gframes[j++ - 1];
					--size;
				}
				gframes[i] = gframes[j] + offset;
			}
			gframes.resize(size);
		};

		if (info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::METIS) {
			if (info::mpi::rank == 0) {
				fixframes();
				info::ecf->input.decomposition.metis_options.continuous = 0;
				edgecut = METIS::call(info::ecf->input.decomposition.metis_options,
						gpartition.size(), gframes.data(), gneighbors.data(),
						0, NULL, NULL, info::mpi::size, gpartition.data());
				info::ecf->input.decomposition.metis_options.continuous = 1;
				profiler::checkpoint("metis");
			}
			profiler::synccheckpoint("synchronize");
			Communication::barrier();
		} else {
			if (MPITools::subset->across.rank == 0) {
				fixframes();
				edistribution = Communication::getDistribution<esint>(gpartition.size(), &MPITools::subset->within);

				switch (info::ecf->input.decomposition.parallel_decomposer) {
				case DecompositionConfiguration::ParallelDecomposer::NONE: break;
				case DecompositionConfiguration::ParallelDecomposer::METIS: break;// never accessed
				case DecompositionConfiguration::ParallelDecomposer::PARMETIS:
					edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
							MPITools::subset->within, edistribution.data(), gframes.data(), gneighbors.data(),
							0, NULL, 0, NULL, NULL, gpartition.data());
					profiler::checkpoint("parmetis");
					break;
				case DecompositionConfiguration::ParallelDecomposer::PTSCOTCH:
					edgecut = PTScotch::call(
							MPITools::subset->within, edistribution.data(), gframes.data(), gneighbors.data(),
							0, NULL, 0, NULL, NULL, gpartition.data());
					profiler::checkpoint("ptscotch");
					break;
				case DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE: break; // never accessed
				case DecompositionConfiguration::ParallelDecomposer::VOLUME: break; // never accessed
				}
			}

			profiler::synccheckpoint("synchronize");
			Communication::barrier(&MPITools::subset->across);
		}

		switch (info::ecf->input.decomposition.parallel_decomposer) {
		case DecompositionConfiguration::ParallelDecomposer::NONE: eslog::checkpointln("MESH: RECLUSTERIZED SKIPPED"); break;
		case DecompositionConfiguration::ParallelDecomposer::METIS: eslog::checkpointln("MESH: RECLUSTERIZED BY METIS"); break;
		case DecompositionConfiguration::ParallelDecomposer::PARMETIS: eslog::checkpointln("MESH: RECLUSTERIZED BY PARMETIS"); break;
		case DecompositionConfiguration::ParallelDecomposer::PTSCOTCH: eslog::checkpointln("MESH: RECLUSTERIZED BY PTSCOTCH"); break;
		case DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE: break; // never accessed
		case DecompositionConfiguration::ParallelDecomposer::VOLUME: break; // never accessed
		}

		if (
				info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::PARMETIS &&
				info::ecf->input.decomposition.parmetis_options.refinement) {

			if (MPITools::subset->across.rank == 0) {
				esint prev = 2 * edgecut;
				while (1.01 * edgecut < prev) {
					prev = edgecut;
					edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_RefineKway,
							MPITools::subset->within, edistribution.data(), gframes.data(), gneighbors.data(),
							0, NULL, 0, NULL, NULL, gpartition.data());
				}
			}
			Communication::barrier(&MPITools::subset->across);
			profiler::synccheckpoint("refinement");
			eslog::checkpointln("MESH: CLUSTERS REFINED BY PARMETIS");
		}

		if (info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::METIS) {
			Communication::broadcastUnknownSize(offsets, &MPITools::singleton->across);
			Communication::scatterv(gpartition, partition, offsets, &MPITools::singleton->across);
			Communication::broadcast(&edgecut, 1, type.mpitype, 0, &MPITools::singleton->across);
		} else {
			Communication::broadcastUnknownSize(offsets, &MPITools::subset->across);
			Communication::scatterv(gpartition, partition, offsets, &MPITools::subset->across);
			Communication::broadcast(&edgecut, 1, type.mpitype, 0, &MPITools::subset->across);
		}


		profiler::synccheckpoint("expand");
		eslog::checkpointln("MESH: MPI PROCESSES EXPANDED");
	}

	profiler::syncend("call_parallel_decomposer");
	return edgecut;
}

void computeElementsClusterization(const ElementStore *elements, const NodeStore *nodes, std::vector<esint> &partition)
{
	if (info::mpi::size == 1) {
		return;
	}

	profiler::syncstart("reclusterize");

	partition.clear();
	partition.resize(elements->epointers->datatarray().size(), info::mpi::rank);

	callParallelDecomposer(elements, nodes, partition);
	profiler::syncend("reclusterize");
}

int getStronglyConnectedComponents(const ElementStore *elements, std::vector<int> &component)
{
	profiler::syncstart("get_strongly_connected_components");
	component.clear();
	component.resize(elements->faceNeighbors->structures(), -1);
	int nextID = 0;
	esint ebegin = elements->offset->datatarray().front();
	esint eend = elements->offset->datatarray().back();
	for (size_t e = 0; e < component.size(); ++e) {
		std::vector<esint> stack;
		if (component[e] == -1) {
			stack.push_back(e);
			component[e] = nextID;
			while (stack.size()) {
				esint current = stack.back();
				stack.pop_back();
				auto faces = elements->faceNeighbors->begin() + current;
				for (auto n = faces->begin(); n != faces->end(); ++n) {
					if (*n != -1 && ebegin <= *n && *n <= eend && component[*n - ebegin] == -1) {
						stack.push_back(*n - ebegin);
						component[*n - ebegin] = nextID;
					}
				}
			}
			nextID++;
		}
	}
	profiler::syncend("get_strongly_connected_components");
	return nextID;
}

void computeContinuousClusterization(const ElementStore *elements, const NodeStore *nodes, const std::vector<esint> &dualDist, const std::vector<esint> &dualData, esint coffset, esint csize, const std::vector<int> &component, const std::vector<int> &neighborsWithMe, std::vector<esint> &partition)
{
	profiler::syncstart("make_continued");

	struct __dinfo__ {
		esint domain, proc, elements, fixed;

		bool operator<(const esint &domain) { return this->domain < domain; }
	};

	std::vector<__dinfo__> sBuffer, gBuffer;
	std::vector<std::vector<__dinfo__> > rBuffer(info::mesh->neighborsWithMe.size());

	sBuffer.resize(csize);
	for (size_t i = 0; i < sBuffer.size(); ++i) {
		sBuffer[i].domain = coffset + i;
		sBuffer[i].proc = info::mpi::rank;
		sBuffer[i].elements = 0;
		sBuffer[i].fixed = i == 0 ? 1 : 0;
	}
	for (size_t i = 0; i < component.size(); ++i) {
		++sBuffer[component[i]].elements;
	}
	for (size_t i = 1, max = 0; i < sBuffer.size(); ++i) {
		if (sBuffer[max].elements < sBuffer[i].elements) {
			sBuffer[max].fixed = 0;
			sBuffer[i].fixed = 1;
			max = i;
		}
	}

	int freeDomain = 1;
	while (freeDomain) {
		freeDomain = 0;
		if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighborsWithMe)) {
			eslog::internalFailure("Cannot exchange domains sizes\n");
		}
		gBuffer.clear();
		for (size_t n = 0; n < rBuffer.size(); ++n) {
			gBuffer.insert(gBuffer.end(), rBuffer[n].begin(), rBuffer[n].end());
		}

		for (esint c = 0; c != csize; ++c) {
			if (!sBuffer[c].fixed) {
				for (esint nc = dualDist[c]; nc != dualDist[c + 1]; ++nc) {
					auto ddesc = std::lower_bound(gBuffer.begin(), gBuffer.end(), dualData[nc]);
					if (ddesc == gBuffer.end() || ddesc->domain != dualData[nc]) {
						eslog::internalFailure("Unknown neighboring domain.\n");
					}
					if (ddesc->fixed) {
						sBuffer[c].proc = ddesc->proc;
						sBuffer[c].fixed = 1;
					} else {
						++freeDomain;
					}
				}
			}
		}
		Communication::allReduce(&freeDomain, NULL, 1, MPI_INT, MPI_SUM);
	}

	profiler::synccheckpoint("decomposition_fixed");
	eslog::checkpointln("MESH: DECOMPOSITION FIXED");

	partition.clear();
	partition.resize(component.size());
	for (size_t i = 0; i < component.size(); ++i) {
		partition[i] = sBuffer[component[i]].proc;
	}

	eslog::startln("MAKE CONTINUED: STARTED", "MAKE CONTINUED");

	eslog::endln("MAKE CONTINUED: FINISHED");
	profiler::syncend("make_continued");
}

void exchangePureElements(ElementStore* &elements, NodeStore* &nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe, const std::vector<esint> &partition)
{
	profiler::syncstart("exchange_pure_elements");
	eslog::startln("EXCHANGE EL: STARTED", "EXCHANGE EL");

	// 0. Compute targets
	// 1. Serialize element data
	// 2. Serialize node data
	// 3. Send serialized data to target (new MPI processes)
	// 4. Deserialize data
	// 5. Balance nodes data to threads
	// 6. Re-index elements (IDs have to be always increasing)

	size_t threads = info::env::OMP_NUM_THREADS;

	// Step 0: Compute targets

	std::vector<int> targets;
	{ // get target processes
		std::vector<std::vector<int> > ttargets(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<int> ttt;
			std::vector<bool> tflags(info::mpi::size, false);
			for (size_t e = elements->offset->datatarray().distribution()[t]; e < elements->offset->datatarray().distribution()[t + 1]; ++e) {
				if (partition[e] != info::mpi::rank) {
					tflags[partition[e]] = true;
				}
			}
			for (int r = 0; r < info::mpi::size; r++) {
				if (tflags[r]) {
					ttt.push_back(r);
				}
			}

			ttargets[t].swap(ttt);
		}

		utils::mergeThreadedUniqueData(ttargets);
		targets = ttargets[0];
	}

	profiler::synccheckpoint("compute_targets");
	eslog::checkpointln("EXCHANGE EL: COMPUTE TARGETS");

	auto t2i = [ & ] (size_t target) {
		return std::lower_bound(targets.begin(), targets.end(), target) - targets.begin();
	};

	ElementStore *newElements = new ElementStore();

	std::vector<std::vector<esint> >    elemsOffset(threads);
	std::vector<std::vector<Element*> > elemsEpointer(threads);
	std::vector<std::vector<esint> >    elemsNodesDistribution(threads);
	std::vector<std::vector<esint> >    elemsNodesData(threads);
	std::vector<std::vector<esint> >    elemsOutputDistribution(threads);
	std::vector<std::vector<esint> >    elemsOutputData(threads);
	std::vector<std::vector<esint> >    elemsRegions(threads);

	NodeStore *newNodes = new NodeStore();

	std::vector<std::vector<esint> >  nodesIDs(threads);
	std::vector<std::vector<Point> >  nodesCoordinates(threads);
	std::vector<std::vector<esint> >  nodesElemsDistribution(threads);
	std::vector<std::vector<esint> >  nodesElemsData(threads);
	std::vector<std::vector<esint> >  nodesOutputDistribution(threads);
	std::vector<std::vector<esint> >  nodesOutputData(threads);
	std::vector<std::vector<esint> >  nodesRegions(threads);

	std::vector<std::vector<std::vector<esint> > >    boundaryEDistribution(boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<esint> > >    boundaryEData(boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<Element*> > > boundaryEPointers(boundaryRegions.size(), std::vector<std::vector<Element*> >(threads));

	// regions are transfered via mask
	int eregionsBitMaskSize = elementsRegions.size() / (8 * sizeof(esint)) + (elementsRegions.size() % (8 * sizeof(esint)) ? 1 : 0);
	int bregionsBitMaskSize = boundaryRegions.size() / (8 * sizeof(esint)) + (boundaryRegions.size() % (8 * sizeof(esint)) ? 1 : 0);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, code, dualsize, dualdata, nodesize, nodeindices, outputOffsetSize, OutputOffsetData)
	std::vector<std::vector<std::vector<esint> > > sElements(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rElements;

	// threads x target x nodes(id, point, linksize, links, regionMask) + size
	std::vector<std::vector<std::vector<esint> > > sNodes(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rNodes;

	// threads x target x boundary(prefix, (code, nodes))
	std::vector<std::vector<std::vector<esint> > > sBoundary(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>(boundaryRegions.size())));
	std::vector<std::vector<esint> > rBoundary;

	// Step 1: Serialize element data

	std::vector<esint> regionElementMask(elements->offset->datatarray().size() * eregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));
			auto begin = std::lower_bound(elementsRegions[r]->elements->datatarray().begin(), elementsRegions[r]->elements->datatarray().end(), elements->offset->datatarray().distribution()[t]);
			auto end = std::lower_bound(elementsRegions[r]->elements->datatarray().begin(), elementsRegions[r]->elements->datatarray().end(), elements->offset->datatarray().distribution()[t + 1]);
			for (auto i = begin; i != end; ++i) {
				regionElementMask[*i * eregionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto offset = elements->offset->datatarray().data();
		auto epointer = elements->epointers->datatarray().data();
		auto enodes = elements->nodes->cbegin(t);
		auto output = elements->inputOffset->cbegin(t);
		auto nIDs = nodes->IDs->datatarray().data();

		std::vector<std::vector<esint> > tsElements(targets.size(), std::vector<esint>({ 0 }));

		std::vector<esint>  telemsOffset;
		std::vector<Element*> telemsEpointer;
		std::vector<esint>  telemsNodesDistribution;
		std::vector<esint>  telemsNodesData;
		std::vector<esint>  telemsNeighborsDistribution;
		std::vector<esint>  telemsNeighborsData;
		std::vector<esint>  telemsOutputDistribution;
		std::vector<esint>  telemsOutputData;
		std::vector<esint>  telemsRegions;
		if (t == 0) {
			telemsNodesDistribution.push_back(0);
			telemsNeighborsDistribution.push_back(0);
			telemsOutputDistribution.push_back(0);
		}

		// estimation
		telemsOffset.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsEpointer.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsNodesDistribution.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsNeighborsDistribution.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsOutputDistribution.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsRegions.reserve(1.5 * elements->offset->datatarray().size() / threads);

		size_t target;
		for (size_t e = elements->offset->datatarray().distribution()[t]; e < elements->offset->datatarray().distribution()[t + 1]; ++e, ++enodes, ++output) {
			PolyElement poly(Element::decode(epointer[e]->code, enodes->size()), enodes->begin());
			if (partition[e] == info::mpi::rank) {
				telemsOffset.push_back(offset[e]);
				telemsEpointer.push_back(epointer[e]);
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					if (poly.isNode(n - enodes->begin())) {
						telemsNodesData.push_back(nIDs[*n]);
					} else {
						telemsNodesData.push_back(*n);
					}
				}
				telemsNodesDistribution.push_back(telemsNodesData.size());
				telemsOutputData.insert(telemsOutputData.end(), output->begin(), output->end());
				telemsOutputDistribution.push_back(telemsOutputData.size());
				telemsRegions.insert(telemsRegions.end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			} else {
				target = t2i(partition[e]);
				tsElements[target].insert(tsElements[target].end(), { offset[e], static_cast<int>(epointer[e]->code) });
				tsElements[target].push_back(enodes->size());
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					if (poly.isNode(n - enodes->begin())) {
						tsElements[target].push_back(nIDs[*n]);
					} else {
						tsElements[target].push_back(*n);
					}
				}
				tsElements[target].push_back(output->size());
				tsElements[target].insert(tsElements[target].end(), output->begin(), output->end());
				tsElements[target].insert(tsElements[target].end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			}
		}

		elemsOffset[t].swap(telemsOffset);
		elemsEpointer[t].swap(telemsEpointer);
		elemsNodesDistribution[t].swap(telemsNodesDistribution);
		elemsNodesData[t].swap(telemsNodesData);
		elemsOutputDistribution[t].swap(telemsOutputDistribution);
		elemsOutputData[t].swap(telemsOutputData);
		elemsRegions[t].swap(telemsRegions);

		sElements[t].swap(tsElements);
	}

	profiler::synccheckpoint("serialize_elements");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE ELEMENTS");

	// Step 2: Serialize node data

	std::vector<esint> regionNodeMask(nodes->size * bregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->nodes) {
				maskOffset = r / (8 * sizeof(esint));
				esint bit = (esint)1 << (r % (8 * sizeof(esint)));
				auto begin = std::lower_bound(boundaryRegions[r]->nodes->datatarray().begin(), boundaryRegions[r]->nodes->datatarray().end(), nodes->IDs->datatarray().distribution()[t]);
				auto end = std::lower_bound(boundaryRegions[r]->nodes->datatarray().begin(), boundaryRegions[r]->nodes->datatarray().end(), nodes->IDs->datatarray().distribution()[t + 1]);
				for (auto i = begin; i != end; ++i) {
					regionNodeMask[*i * bregionsBitMaskSize + maskOffset] |= bit;
				}
			}
		}
	}

	esint eBegin = elements->offset->datatarray().size();
	Communication::exscan(eBegin);
	esint eEnd = eBegin + elements->offset->datatarray().size();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = nodes->IDs->datatarray();
		const auto &coordinates = nodes->coordinates->datatarray();
		auto elems = nodes->elements->cbegin(t);
		auto output = nodes->inputOffset->cbegin(t);

		std::vector<esint>  tnodesIDs;
		std::vector<Point>  tnodesCoordinates;
		std::vector<esint>  tnodesElemsDistribution;
		std::vector<esint>  tnodesElemsData;
		std::vector<esint>  tnodesOutputDistribution;
		std::vector<esint>  tnodesOutputData;
		std::vector<esint>  tnodesRegions;

		if (t == 0) {
			tnodesElemsDistribution.push_back(0);
			tnodesOutputDistribution.push_back(0);
		}

		std::vector<std::vector<esint> > tsNodes(targets.size(), std::vector<esint>({ 0 }));

		tnodesIDs.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesCoordinates.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesElemsDistribution.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesOutputDistribution.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesRegions.reserve(1.5 * nodes->IDs->datatarray().size() / threads);

		size_t target;
		std::vector<bool> last(targets.size() + 1); // targets + me
		for (size_t n = nodes->IDs->datatarray().distribution()[t]; n < nodes->IDs->datatarray().distribution()[t + 1]; ++n, ++elems, ++output) {
			std::fill(last.begin(), last.end(), false);
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				if (eBegin <= *e && *e < eEnd) {
					target = t2i(partition[*e - eBegin]);
					if (!last[target] && partition[*e - eBegin] != info::mpi::rank) {
						tsNodes[target].push_back(IDs[n]);
						tsNodes[target].insert(tsNodes[target].end(), reinterpret_cast<const esint*>(coordinates.data() + n), reinterpret_cast<const esint*>(coordinates.data() + n + 1));
						tsNodes[target].push_back(elems->size());
						tsNodes[target].insert(tsNodes[target].end(), elems->begin(), elems->end());
						tsNodes[target].push_back(output->size());
						tsNodes[target].insert(tsNodes[target].end(), output->begin(), output->end());
						tsNodes[target].insert(tsNodes[target].end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last[target] = true;
					}
					if (!last.back() && partition[*e - eBegin] == info::mpi::rank) {
						tnodesIDs.push_back(IDs[n]);
						tnodesCoordinates.push_back(coordinates[n]);
						tnodesElemsData.insert(tnodesElemsData.end(), elems->begin(), elems->end());
						tnodesElemsDistribution.push_back(tnodesElemsData.size());
						tnodesOutputData.insert(tnodesOutputData.end(), output->begin(), output->end());
						tnodesOutputDistribution.push_back(tnodesOutputData.size());
						tnodesRegions.insert(tnodesRegions.end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last.back() = true;
					}
				}
			}
		}

		nodesIDs[t].swap(tnodesIDs);
		nodesCoordinates[t].swap(tnodesCoordinates);
		nodesElemsDistribution[t].swap(tnodesElemsDistribution);
		nodesElemsData[t].swap(tnodesElemsData);
		nodesOutputDistribution[t].swap(tnodesOutputDistribution);
		nodesOutputData[t].swap(tnodesOutputData);
		nodesRegions[t].swap(tnodesRegions);

		sNodes[t].swap(tsNodes);
	}

	profiler::synccheckpoint("serialize_nodes");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE NODES");

	// Step 2.1: Serialize boundary regions data

	std::vector<esint> emembership;

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->originalDimension) {
			emembership.clear();
			emembership.resize(boundaryRegions[r]->epointers->datatarray().size());
			const std::vector<size_t> &distribution = boundaryRegions[r]->epointers->datatarray().distribution();

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto enodes = boundaryRegions[r]->elements->cbegin() + distribution[t];
				auto ep = boundaryRegions[r]->epointers->datatarray().begin(t);
				std::vector<esint> nlinks;
				int counter;
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes, ++ep) {
					PolyElement poly(Element::decode((*ep)->code, enodes->size()), enodes->begin());
					nlinks.clear();
					for (auto n = enodes->begin(); n != enodes->end(); ++n) {
						if (poly.isNode(n - enodes->begin())) {
							auto links = nodes->elements->cbegin() + *n;
							nlinks.insert(nlinks.end(), links->begin(), links->end());
						}
					}
					std::sort(nlinks.begin(), nlinks.end());
					counter = 1;
					for (size_t i = 1; i < nlinks.size(); ++i) {
						if (nlinks[i - 1] == nlinks[i]) {
							++counter;
							if (counter == poly.size && eBegin <= nlinks[i]) {
								emembership[e] = nlinks[i] - eBegin;
								break;
							}
						} else {
							counter = 1;
						}
					}
				}
			}

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				std::vector<esint>  tboundaryEDistribution;
				std::vector<esint>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;
				if (t == 0) {
					tboundaryEDistribution.push_back(0);
				}

				std::vector<size_t> tsize(targets.size());
				for (size_t i = 0; i < targets.size(); i++) {
					tsize[i] = sBoundary[t][i].size();
				}
				esint target;
				auto enodes = boundaryRegions[r]->elements->cbegin() + distribution[t];
				const auto &IDs = nodes->IDs->datatarray();
				const auto &epointer = boundaryRegions[r]->epointers->datatarray();
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					PolyElement poly(Element::decode(epointer[e]->code, enodes->size()), enodes->begin());
					if (partition[emembership[e]] == info::mpi::rank) {
						tboundaryEPointers.push_back(epointer[e]);
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							if (poly.isNode(n - enodes->begin())) {
								tboundaryEData.push_back(IDs[*n]);
							} else {
								tboundaryEData.push_back(*n);
							}
						}
						tboundaryEDistribution.push_back(tboundaryEData.size());
					} else {
						target = t2i(partition[emembership[e]]);
						sBoundary[t][target].push_back(static_cast<int>(epointer[e]->code));
						sBoundary[t][target].push_back(enodes->size());
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							if (poly.isNode(n - enodes->begin())) {
								sBoundary[t][target].push_back(IDs[*n]);
							} else {
								sBoundary[t][target].push_back(*n);
							}
						}
					}
				}
				for (size_t i = 0; i < targets.size(); i++) {
					sBoundary[t][i][r] = sBoundary[t][i].size() - tsize[i];
				}

				boundaryEDistribution[r][t].swap(tboundaryEDistribution);
				boundaryEData[r][t].swap(tboundaryEData);
				boundaryEPointers[r][t].swap(tboundaryEPointers);
			}
		}
	}

	profiler::synccheckpoint("serialize_boundary");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE BOUNDARIES");

	// Step 3: Send data to target processes

	#pragma omp parallel for
	for (size_t target = 0; target < targets.size(); ++target) {
		sElements[0][target].front() = sElements[0][target].size();
		sNodes[0][target].front() = sNodes[0][target].size();
		for (size_t t = 1; t < threads; t++) {
			sElements[t][target].front() = sElements[t][target].size();
			sElements[0][target].insert(sElements[0][target].end(), sElements[t][target].begin(), sElements[t][target].end());
			sNodes[t][target].front() = sNodes[t][target].size();
			sNodes[0][target].insert(sNodes[0][target].end(), sNodes[t][target].begin(), sNodes[t][target].end());
			sBoundary[0][target].insert(sBoundary[0][target].end(), sBoundary[t][target].begin(), sBoundary[t][target].end());
		}
	}

	if (!Communication::sendVariousTargets(sElements[0], rElements, targets)) {
		eslog::internalFailure("exchange elements data.\n");
	}

	if (!Communication::sendVariousTargets(sNodes[0], rNodes, targets)) {
		eslog::internalFailure("exchange nodes data.\n");
	}

	if (!Communication::sendVariousTargets(sBoundary[0], rBoundary, targets)) {
		eslog::internalFailure("exchange boundary data.\n");
	}

	profiler::synccheckpoint("exchange");
	eslog::checkpointln("EXCHANGE EL: EXCHANGE");

	// Step 4: Deserialize element data
	for (size_t i = 0; i < rElements.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rElements[i].size()) {
			rdistribution.push_back(p += rElements[i][p]);
		}

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<esint>  telemsOffset;
			std::vector<Element*> telemsEpointer;
			std::vector<esint>  telemsNodesDistribution;
			std::vector<esint>  telemsNodesData;
			std::vector<esint>  telemsOutputDistribution;
			std::vector<esint>  telemsOutputData;
			std::vector<esint>  telemsRegions;

			telemsOffset.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsEpointer.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsNodesDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsOutputDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsRegions.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0, outputOffset = 0;

			if (elemsNodesDistribution[t].size()) {
				distOffset = elemsNodesDistribution[t].back();
			}
			if (elemsOutputDistribution[t].size()) {
				outputOffset = elemsOutputDistribution[t].back();
			}
			if (t == 0 && elemsNodesDistribution[t].size() == 0) {
				telemsNodesDistribution.push_back(0);
			}
			if (t == 0 && elemsOutputDistribution[t].size() == 0) {
				telemsOutputDistribution.push_back(0);
			}

			for (size_t e = rdistribution[t] + 1; e < rdistribution[t + 1]; ) {
				telemsOffset.push_back(rElements[i][e++]);
				telemsEpointer.push_back(&Mesh::edata[rElements[i][e++]]);
				telemsNodesData.insert(telemsNodesData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNodesDistribution.push_back(telemsNodesData.size() + distOffset);
				e += rElements[i][e++]; // nodes + nodes size
				telemsOutputData.insert(telemsOutputData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsOutputDistribution.push_back(telemsOutputData.size() + outputOffset);
				e += rElements[i][e++]; // output + output size
				telemsRegions.insert(telemsRegions.end(), rElements[i].begin() + e, rElements[i].begin() + e + eregionsBitMaskSize);
				e += eregionsBitMaskSize;
			}

			elemsOffset[t].insert(elemsOffset[t].end(), telemsOffset.begin(), telemsOffset.end());
			elemsEpointer[t].insert(elemsEpointer[t].end(), telemsEpointer.begin(), telemsEpointer.end());
			elemsNodesDistribution[t].insert(elemsNodesDistribution[t].end(), telemsNodesDistribution.begin(), telemsNodesDistribution.end());
			elemsNodesData[t].insert(elemsNodesData[t].end(), telemsNodesData.begin(), telemsNodesData.end());
			elemsOutputDistribution[t].insert(elemsOutputDistribution[t].end(), telemsOutputDistribution.begin(), telemsOutputDistribution.end());
			elemsOutputData[t].insert(elemsOutputData[t].end(), telemsOutputData.begin(), telemsOutputData.end());
			elemsRegions[t].insert(elemsRegions[t].end(), telemsRegions.begin(), telemsRegions.end());
		}
	}

	profiler::synccheckpoint("deserialize_elements");
	eslog::checkpointln("EXCHANGE EL: DESERIALIZE ELEMENTS");

	// Step 4: Deserialize node data
	std::vector<esint> nodeset;
	for (size_t t = 0; t < threads; t++) {
		nodeset.insert(nodeset.end(), nodesIDs[t].begin(), nodesIDs[t].end());
	}
	std::sort(nodeset.begin(), nodeset.end());
	for (size_t i = 0; i < rNodes.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rNodes[i].size()) {
			rdistribution.push_back(p += rNodes[i][p]);
		}
		std::vector<std::vector<esint> > tnodeset(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point point;

			std::vector<esint>  tnodesIDs;
			std::vector<Point>  tnodesCoordinates;
			std::vector<esint>  tnodesElemsDistribution;
			std::vector<esint>  tnodesElemsData;
			std::vector<esint>  tnodesOutputDistribution;
			std::vector<esint>  tnodesOutputData;
			std::vector<esint>  tnodesRegions;
			std::vector<esint>  tnodeSet;
			std::vector<esint>  tnpermutation;

			tnodesIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesCoordinates.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesElemsDistribution.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesOutputDistribution.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesRegions.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodeSet.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0, outputOffset = 0;

			if (nodesElemsDistribution[t].size()) {
				distOffset = nodesElemsDistribution[t].back();
			}
			if (nodesOutputDistribution[t].size()) {
				outputOffset = nodesOutputDistribution[t].back();
			}
			if (t == 0 && nodesElemsDistribution[t].size() == 0) {
				tnodesElemsDistribution.push_back(0);
			}
			if (t == 0 && nodesOutputDistribution[t].size() == 0) {
				tnodesOutputDistribution.push_back(0);
			}

			for (size_t n = rdistribution[t] + 1; n < rdistribution[t + 1]; ) {
				tnpermutation.push_back(n);
				n += 1 + sizeof(Point) / sizeof(esint); // id, Point
				n += 1 + rNodes[i][n]; // linksize, links
				n += 1 + rNodes[i][n]; // outputsize, outputs
				n += bregionsBitMaskSize; // region mask
			}
			std::sort(tnpermutation.begin(), tnpermutation.end(), [&] (esint n1, esint n2) {
				return rNodes[i][n1] < rNodes[i][n2];
			});
			auto nodesetit = nodeset.begin();
			if (rdistribution[t] + 1 < rdistribution[t + 1]) {
				nodesetit = std::lower_bound(nodeset.begin(), nodeset.end(), rNodes[i][tnpermutation.front()]);
			}
			size_t index;
			for (size_t n = 0; n < tnpermutation.size(); n++) {
				index = tnpermutation[n];
				while (nodesetit != nodeset.end() && *nodesetit < rNodes[i][index]) ++nodesetit;
				if (nodesetit == nodeset.end() || *nodesetit != rNodes[i][index]) {
					tnodesIDs.push_back(rNodes[i][index]);
					tnodeSet.push_back(rNodes[i][index]);
					index += 1; //ID
					memcpy(reinterpret_cast<void*>(&point), rNodes[i].data() + index, sizeof(Point));
					tnodesCoordinates.push_back(point);
					index += sizeof(Point) / sizeof(esint); // points

					tnodesElemsData.insert(tnodesElemsData.end(), rNodes[i].begin() + index + 1, rNodes[i].begin() + index + 1 + rNodes[i][index]);
					tnodesElemsDistribution.push_back(tnodesElemsData.size() + distOffset);
					index += rNodes[i][index] + 1; // linksize + links

					tnodesOutputData.insert(tnodesOutputData.end(), rNodes[i].begin() + index + 1, rNodes[i].begin() + index + 1 + rNodes[i][index]);
					tnodesOutputDistribution.push_back(tnodesOutputData.size() + outputOffset);
					index += rNodes[i][index] + 1; // outputsize + output

					tnodesRegions.insert(tnodesRegions.end(), rNodes[i].begin() + index, rNodes[i].begin() + index + bregionsBitMaskSize);
					index += bregionsBitMaskSize; // region mask
				}
			}

			nodesIDs[t].insert(nodesIDs[t].end(), tnodesIDs.begin(), tnodesIDs.end());
			nodesCoordinates[t].insert(nodesCoordinates[t].end(), tnodesCoordinates.begin(), tnodesCoordinates.end());
			nodesElemsDistribution[t].insert(nodesElemsDistribution[t].end(), tnodesElemsDistribution.begin(), tnodesElemsDistribution.end());
			nodesElemsData[t].insert(nodesElemsData[t].end(), tnodesElemsData.begin(), tnodesElemsData.end());
			nodesOutputDistribution[t].insert(nodesOutputDistribution[t].end(), tnodesOutputDistribution.begin(), tnodesOutputDistribution.end());
			nodesOutputData[t].insert(nodesOutputData[t].end(), tnodesOutputData.begin(), tnodesOutputData.end());
			nodesRegions[t].insert(nodesRegions[t].end(), tnodesRegions.begin(), tnodesRegions.end());
			tnodeset[t].swap(tnodeSet);
		}

		size_t nsize = nodeset.size();
		for (size_t t = 0; t < threads; t++) {
			nodeset.insert(nodeset.end(), tnodeset[t].begin(), tnodeset[t].end());
		}
		std::inplace_merge(nodeset.begin(), nodeset.begin() + nsize, nodeset.end());
	}

	utils::threadDistributionToFullDistribution(elemsNodesDistribution);
	utils::threadDistributionToFullDistribution(elemsOutputDistribution);
	utils::threadDistributionToFullDistribution(nodesElemsDistribution);
	utils::threadDistributionToFullDistribution(nodesOutputDistribution);

	profiler::synccheckpoint("deserialize_nodes");
	eslog::checkpointln("EXCHANGE EL: DESERIALIZE NODES");

	// Step 4: Deserialize boundary data
	for (size_t n = 0; n < rBoundary.size(); ++n) {
		std::vector<std::vector<esint> > toffset(boundaryRegions.size()), tsize(boundaryRegions.size());
		esint offset = 0, p = 0;
		for (size_t t = 0; t < threads; t++) {
			for (size_t r = 0; r < boundaryRegions.size(); r++) {
				toffset[r].push_back(offset + boundaryRegions.size());
				tsize[r].push_back(rBoundary[n][p + r]);
				offset += rBoundary[n][p + r];
			}
			offset += boundaryRegions.size();
			p = offset;
		}
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {

				std::vector<esint>  tboundaryEDistribution;
				std::vector<esint>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;

				esint distOffset = 0;

				if (boundaryEDistribution[r][t].size()) {
					distOffset = boundaryEDistribution[r][t].back();
				}
				if (t == 0 && boundaryEDistribution[r][t].size() == 0) {
					tboundaryEDistribution.push_back(0);
				}

				for (esint i = toffset[r][t]; i < toffset[r][t] + tsize[r][t];) {
					tboundaryEPointers.push_back(&Mesh::edata[rBoundary[n][i++]]);
					tboundaryEData.insert(tboundaryEData.end(), rBoundary[n].begin() + i + 1, rBoundary[n].begin() + i + 1 + rBoundary[n][i]);
					tboundaryEDistribution.push_back(tboundaryEData.size() + distOffset);
					i += rBoundary[n][i] + 1;
				}

				boundaryEPointers[r][t].insert(boundaryEPointers[r][t].end(), tboundaryEPointers.begin(), tboundaryEPointers.end());
				boundaryEDistribution[r][t].insert(boundaryEDistribution[r][t].end(), tboundaryEDistribution.begin(), tboundaryEDistribution.end());
				boundaryEData[r][t].insert(boundaryEData[r][t].end(), tboundaryEData.begin(), tboundaryEData.end());
			}
		}
	}

	profiler::synccheckpoint("deserialize_boundary");
	eslog::checkpointln("EXCHANGE EL: DESERIALIZE BOUDARIES");

	// elements are redistributed later while decomposition -> distribution is not changed now
	std::vector<size_t> elemDistribution(threads + 1);
	for (size_t t = 1; t <= threads; t++) {
		elemDistribution[t] = elemDistribution[t - 1] + elemsOffset[t - 1].size();
	}

	newElements->offset = new serializededata<esint, esint>(1, elemsOffset);
	newElements->epointers = new serializededata<esint, Element*>(1, elemsEpointer);
	newElements->nodes = new serializededata<esint, esint>(elemsNodesDistribution, elemsNodesData); // global IDs

	newElements->regions = new serializededata<esint, esint>(eregionsBitMaskSize, elemsRegions);
	newElements->inputOffset = new serializededata<esint, esint>(elemsOutputDistribution, elemsOutputData);

	newElements->distribution.process.size = newElements->offset->structures();
	newElements->distribution.threads = newElements->offset->datatarray().distribution();

	// Step 5: Balance node data to threads
	std::vector<size_t> nodeDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		nodeDistribution[t] = nodeDistribution[t - 1] + nodesIDs[t - 1].size();
	}

	serializededata<esint, esint>::balance(1, nodesIDs);
	serializededata<esint, Point>::balance(1, nodesCoordinates);
	serializededata<esint, esint>::balance(nodesElemsDistribution, nodesElemsData);

	newNodes->IDs = new serializededata<esint, esint>(1, nodesIDs);
	newNodes->coordinates = new serializededata<esint, Point>(1, nodesCoordinates);
	newNodes->elements = new serializededata<esint, esint>(nodesElemsDistribution, nodesElemsData);
	newNodes->inputOffset = new serializededata<esint, esint>(nodesOutputDistribution, nodesOutputData);
	newNodes->size = newNodes->IDs->datatarray().size();
	newNodes->distribution = newNodes->IDs->datatarray().distribution();

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->originalDimension) {
			delete boundaryRegions[r]->elements;
			delete boundaryRegions[r]->epointers;

			utils::threadDistributionToFullDistribution(boundaryEDistribution[r]);
			boundaryRegions[r]->elements = new serializededata<esint, esint>(boundaryEDistribution[r], boundaryEData[r]);
			boundaryRegions[r]->epointers = new serializededata<esint, Element*>(1, boundaryEPointers[r]);
		}
	}

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		esint maskOffset = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));
		delete elementsRegions[r]->elements;
		std::vector<std::vector<esint> > regionelems(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = 0; i < elemsRegions[t].size(); i += eregionsBitMaskSize) {
				if (elemsRegions[t][i + maskOffset] & bit) {
					regionelems[t].push_back(elemDistribution[t] + i / eregionsBitMaskSize);
				}
			}
		}

		serializededata<esint, esint>::balance(1, regionelems);
		elementsRegions[r]->elements = new serializededata<esint, esint>(1, regionelems);
	}

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->nodes) {
			esint maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));
			delete boundaryRegions[r]->nodes;
			std::vector<std::vector<esint> > regionnodes(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = 0; i < nodesRegions[t].size(); i += bregionsBitMaskSize) {
					if (nodesRegions[t][i + maskOffset] & bit) {
						regionnodes[t].push_back(nodeDistribution[t] + i / bregionsBitMaskSize);
					}
				}
			}

			serializededata<esint, esint>::balance(1, regionnodes);
			boundaryRegions[r]->nodes = new serializededata<esint, esint>(1, regionnodes);
		}
	}

	std::vector<size_t> eIDsOLD = Communication::getDistribution(elements->epointers->datatarray().size());
	std::vector<size_t> eIDsNEW = Communication::getDistribution(newElements->epointers->datatarray().size());

	for (size_t t = 1; t < threads; ++t) {
		elemsOffset[0].insert(elemsOffset[0].end(), elemsOffset[t].begin(), elemsOffset[t].end());
	}

	std::vector<esint> epermutation(elemsOffset[0].size());
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return elemsOffset[0][i] < elemsOffset[0][j]; });

	std::vector<esint> sortedElements(epermutation.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		for (size_t e = elemDistribution[t]; e < elemDistribution[t + 1]; e++) {
			sortedElements[e] = elemsOffset[0][epermutation[e]];
		}
		utils::sortAndRemoveDuplicates(nodesElemsData[t]);
	}
	utils::inplaceMerge(nodesElemsData);
	utils::removeDuplicates(nodesElemsData[0]);

	std::vector<std::vector<esint> > requestedIDs, receivedTargets, IDrequests;
	std::vector<std::vector<int> > IDtargets(threads);
	std::vector<int> sources;

	std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, info::mpi::size);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		std::vector<int> ttargets;

		auto eIDbegin = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[rdistribution[t]]);
		auto eIDend   = std::lower_bound(eIDbegin, nodesElemsData[0].end(), eIDsOLD[rdistribution[t + 1]]);

		for (size_t r = rdistribution[t]; r < rdistribution[t + 1]; r++) {
			eIDbegin = std::lower_bound(eIDbegin, eIDend, eIDsOLD[r]);
			auto end = std::lower_bound(eIDbegin, eIDend, eIDsOLD[r + 1]);
			if (eIDbegin != end) {
				ttargets.push_back(r);
			}
			eIDbegin = end;
		}

		IDtargets[t].swap(ttargets);
	}

	for (size_t t = 1; t < threads; ++t) {
		IDtargets[0].insert(IDtargets[0].end(), IDtargets[t].begin(), IDtargets[t].end());
	}

	requestedIDs.resize(IDtargets[0].size());

	#pragma omp parallel for
	for (size_t t = 0; t < IDtargets[0].size(); t++) {
		auto mybegin = std::lower_bound(sortedElements.begin(), sortedElements.end(), eIDsOLD[IDtargets[0][t]]);
		auto myend   = std::lower_bound(sortedElements.begin(), sortedElements.end(), eIDsOLD[IDtargets[0][t] + 1]);
		auto nbegin = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[IDtargets[0][t]]);
		auto nend   = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[IDtargets[0][t] + 1]);
		requestedIDs[t].resize(nend - nbegin);
		requestedIDs[t].resize(std::set_difference(nbegin, nend, mybegin, myend, requestedIDs[t].begin()) - requestedIDs[t].begin());
	}

	profiler::synccheckpoint("postprocessing");
	eslog::checkpointln("EXCHANGE EL: POST-PROCESSING");

	if (!Communication::sendVariousTargets(requestedIDs, IDrequests, IDtargets[0], sources)) {
		eslog::internalFailure("exchange ID requests.\n");
	}

	for (size_t r = 0; r < IDrequests.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, IDrequests[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			if (rdistribution[t] != rdistribution[t + 1]) {
				for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
					IDrequests[r][e] = partition[IDrequests[r][e] - eIDsOLD[info::mpi::rank]];
				}
			}
		}
	}

	if (!Communication::sendVariousTargets(IDrequests, receivedTargets, sources)) {
		eslog::internalFailure("return ID targets.\n");
	}

	profiler::synccheckpoint("idrequests");
	eslog::checkpointln("EXCHANGE EL: IDS REQUESTS");

	IDtargets.clear();
	IDtargets.resize(threads);

	std::vector<std::vector<bool> > newtargets(threads, std::vector<bool>(info::mpi::size, false));
	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				newtargets[t][receivedTargets[r][e]] = true;
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		std::vector<int> ttargets;
		for (int r = 0; r < info::mpi::size; r++) {
			if (r != info::mpi::rank && newtargets[t][r]) {
				ttargets.push_back(r);
			}
		}
		IDtargets[t].swap(ttargets);
	}
	for (size_t t = 1; t < threads; ++t) {
		IDtargets[0].insert(IDtargets[0].end(), IDtargets[t].begin(), IDtargets[t].end());
	}
	utils::sortAndRemoveDuplicates(IDtargets[0]);

	std::vector<std::vector<std::vector<esint> > > newIDrequests(threads, std::vector<std::vector<esint> >(IDtargets[0].size()));
	std::vector<std::vector<esint> > newIDs;

	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			size_t tindex;
			std::vector<std::vector<esint> > tnewIDrequests(IDtargets[0].size());

			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				tindex = std::lower_bound(IDtargets[0].begin(), IDtargets[0].end(), receivedTargets[r][e]) - IDtargets[0].begin();
				tnewIDrequests[tindex].push_back(requestedIDs[r][e]);
			}

			for (size_t n = 0; n < IDtargets[0].size(); n++) {
				newIDrequests[t][n].insert(newIDrequests[t][n].end(), tnewIDrequests[n].begin(), tnewIDrequests[n].end());
			}
		}

		for (size_t t = 1; t < threads; ++t) {
			for (size_t n = 0; n < IDtargets[0].size(); n++) {
				newIDrequests[0][n].insert(newIDrequests[0][n].end(), newIDrequests[t][n].begin(), newIDrequests[t][n].end());
				newIDrequests[t][n].clear();
			}
		}
	}

	profiler::synccheckpoint("process_requests");
	eslog::checkpointln("EXCHANGE EL: PROC REQUESTS");

	if (!Communication::sendVariousTargets(newIDrequests[0], IDrequests, IDtargets[0], sources)) {
		eslog::internalFailure("request new ID.\n");
	}

	for (size_t r = 0; r < IDrequests.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, IDrequests[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			if (rdistribution[t] != rdistribution[t + 1]) {
				auto eit = std::lower_bound(sortedElements.begin(), sortedElements.end(), IDrequests[r][rdistribution[t]]);
				for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
					while (eit != sortedElements.end() && *eit < IDrequests[r][e]) ++eit;
					IDrequests[r][e] = epermutation[eit - sortedElements.begin()] + eIDsNEW[info::mpi::rank];
				}
			}
		}
	}

	if (!Communication::sendVariousTargets(IDrequests, newIDs, sources)) {
		eslog::internalFailure("return new ID.\n");
	}

	profiler::synccheckpoint("exchange_ids");
	eslog::checkpointln("EXCHANGE EL: NEW IDS");

	std::vector<size_t> offsets = { 0 };
	for (size_t r = 0; r < newIDs.size(); r++) {
		offsets.push_back(offsets.back() + newIDs[r].size());
	}
	std::vector<std::pair<esint, esint> > IDMap(offsets.back());

	for (size_t r = 0; r < newIDs.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, newIDrequests[0][r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				IDMap[offsets[r] + e].first = newIDrequests[0][r][e];
				IDMap[offsets[r] + e].second = newIDs[r][e];
			}
		}
	}

	utils::sortWithInplaceMerge(IDMap, offsets);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto elem = newNodes->elements->begin(t); elem != newNodes->elements->end(t); ++elem) {
			for (auto e = elem->begin(); e != elem->end(); ++e) {
				auto mapit = std::lower_bound(IDMap.begin(), IDMap.end(), std::make_pair(*e, (esint)0));
				if (mapit == IDMap.end() || mapit->first != *e) {
					*e = epermutation[std::lower_bound(sortedElements.begin(), sortedElements.end(), *e) - sortedElements.begin()] + eIDsNEW[info::mpi::rank];
				} else {
					*e = mapit->second;
				}
			}
			std::sort(elem->begin(), elem->end());
		}
	}

	std::vector<esint> permutation(newNodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return newNodes->IDs->datatarray()[i] < newNodes->IDs->datatarray()[j]; });

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ep = newElements->epointers->datatarray().begin(t);
		for (auto enodes = newElements->nodes->begin(t); enodes != newElements->nodes->end(t); ++enodes, ++ep) {
			PolyElement poly(Element::decode((*ep)->code, enodes->size()), enodes->begin());
			for (size_t n = 0; n < enodes->size(); ++n) {
				if (poly.isNode(n)) {
					enodes->at(n) = *std::lower_bound(permutation.begin(), permutation.end(), enodes->at(n), [&] (esint i, esint val) {
						return newNodes->IDs->datatarray()[i] < val;
					});
				}
			}
		}
	}

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->elements) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto ep = boundaryRegions[r]->epointers->datatarray().begin(t);
				for (auto enodes = boundaryRegions[r]->elements->begin(t); enodes != boundaryRegions[r]->elements->end(t); ++enodes, ++ep) {
					PolyElement poly(Element::decode((*ep)->code, enodes->size()), enodes->begin());
					for (size_t n = 0; n < enodes->size(); ++n) {
						if (poly.isNode(n)) {
							enodes->at(n) = *std::lower_bound(permutation.begin(), permutation.end(), enodes->at(n), [&] (esint i, esint val) {
								return newNodes->IDs->datatarray()[i] < val;
							});
						}
					}
				}
			}
		}
	}

	std::vector<std::vector<esint> > rankBoundaries(threads);
	std::vector<std::vector<int> > rankData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> trankBoundaries;
		std::vector<int> trankData;
		if (t == 0) {
			trankBoundaries.push_back(0);
		}

		int rank;
		for (auto elem = newNodes->elements->begin(t); elem != newNodes->elements->end(t); ++elem) {
			trankData.push_back(std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *elem->begin() + 1) - eIDsNEW.begin() - 1);
			for (auto e = elem->begin() + 1; e != elem->end(); ++e) {
				rank = std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *e + 1) - eIDsNEW.begin() - 1;
				if (rank != trankData.back()) {
					trankData.push_back(rank);
				}
			}
			trankBoundaries.push_back(trankData.size());
		}

		rankBoundaries[t].swap(trankBoundaries);
		rankData[t].swap(trankData);
	}

	utils::threadDistributionToFullDistribution(rankBoundaries);

	newNodes->ranks = new serializededata<esint, int>(rankBoundaries, rankData);

	std::iota(newElements->offset->datatarray().begin(), newElements->offset->datatarray().end(), eIDsNEW[info::mpi::rank]);
	newElements->distribution.process.offset = eIDsNEW[info::mpi::rank];
	newElements->distribution.process.totalSize = elements->distribution.process.totalSize;
	std::swap(elements, newElements);
	std::swap(nodes, newNodes);
	neighbors.clear();
	for (size_t t = 0; t < IDtargets[0].size(); t++) {
		if (IDtargets[0][t] != info::mpi::rank) {
			neighbors.push_back(IDtargets[0][t]);
		}
	}
	neighborsWithMe = neighbors;
	neighborsWithMe.push_back(info::mpi::rank);
	std::sort(neighborsWithMe.begin(), neighborsWithMe.end());

	delete newElements;
	delete newNodes;

	profiler::synccheckpoint("finish");
	profiler::syncend("exchange_pure_elements");
	eslog::endln("EXCHANGE EL: FINISH");
	eslog::checkpointln("MESH: ELEMENTS EXCHANGED");
}

void exchangeElements(ElementStore* &elements, NodeStore* &nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe, const std::vector<esint> &partition)
{
	profiler::syncstart("exchange_elements");
	eslog::startln("EXCHANGE EL: STARTED", "EXCHANGE EL");

	// 0. Compute targets
	// 1. Serialize element data
	// 2. Serialize node data
	// 3. Send serialized data to target (new MPI processes)
	// 4. Deserialize data
	// 5. Balance nodes data to threads
	// 6. Re-index elements (IDs have to be always increasing)

	size_t threads = info::env::OMP_NUM_THREADS;

	// Step 0: Compute targets

	std::vector<int> targets;
	{ // get target processes
		std::vector<std::vector<int> > ttargets(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<int> ttt;
			std::vector<bool> tflags(info::mpi::size, false);
			for (size_t e = elements->offset->datatarray().distribution()[t]; e < elements->offset->datatarray().distribution()[t + 1]; ++e) {
				if (partition[e] != info::mpi::rank) {
					tflags[partition[e]] = true;
				}
			}
			for (int r = 0; r < info::mpi::size; r++) {
				if (tflags[r]) {
					ttt.push_back(r);
				}
			}

			ttargets[t].swap(ttt);
		}

		utils::mergeThreadedUniqueData(ttargets);
		targets = ttargets[0];
	}

	profiler::synccheckpoint("compute_targets");
	eslog::checkpointln("EXCHANGE EL: COMPUTE TARGETS");

	auto t2i = [ & ] (size_t target) {
		return std::lower_bound(targets.begin(), targets.end(), target) - targets.begin();
	};

	ElementStore *newElements = new ElementStore();

	std::vector<std::vector<esint> >    elemsOffset(threads);
	std::vector<std::vector<Element*> > elemsEpointer(threads);
	std::vector<std::vector<esint> >    elemsNodesDistribution(threads);
	std::vector<std::vector<esint> >    elemsNodesData(threads);
	std::vector<std::vector<esint> >    elemsNeighborsDistribution(threads);
	std::vector<std::vector<esint> >    elemsNeighborsData(threads);
	std::vector<std::vector<esint> >    elemsOutputDistribution(threads);
	std::vector<std::vector<esint> >    elemsOutputData(threads);
	std::vector<std::vector<esint> >    elemsRegions(threads);

	NodeStore *newNodes = new NodeStore();

	std::vector<std::vector<esint> >  nodesIDs(threads);
	std::vector<std::vector<Point> >  nodesCoordinates(threads);
	std::vector<std::vector<esint> >  nodesElemsDistribution(threads);
	std::vector<std::vector<esint> >  nodesElemsData(threads);
	std::vector<std::vector<esint> >  nodesOutputDistribution(threads);
	std::vector<std::vector<esint> >  nodesOutputData(threads);
	std::vector<std::vector<esint> >  nodesRegions(threads);

	std::vector<std::vector<std::vector<esint> > >    boundaryEDistribution(boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<esint> > >    boundaryEData(boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<Element*> > > boundaryEPointers(boundaryRegions.size(), std::vector<std::vector<Element*> >(threads));

	// regions are transfered via mask
	int eregionsBitMaskSize = elementsRegions.size() / (8 * sizeof(esint)) + (elementsRegions.size() % (8 * sizeof(esint)) ? 1 : 0);
	int bregionsBitMaskSize = boundaryRegions.size() / (8 * sizeof(esint)) + (boundaryRegions.size() % (8 * sizeof(esint)) ? 1 : 0);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, code, dualsize, dualdata, nodesize, nodeindices, outputOffsetSize, OutputOffsetData)
	std::vector<std::vector<std::vector<esint> > > sElements(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rElements;

	// threads x target x nodes(id, point, linksize, links, regionMask) + size
	std::vector<std::vector<std::vector<esint> > > sNodes(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rNodes;

	// threads x target x boundary(prefix, (code, nodes))
	std::vector<std::vector<std::vector<esint> > > sBoundary(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>(boundaryRegions.size())));
	std::vector<std::vector<esint> > rBoundary;

	// Step 1: Serialize element data

	std::vector<esint> regionElementMask(elements->offset->datatarray().size() * eregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));
			auto begin = std::lower_bound(elementsRegions[r]->elements->datatarray().begin(), elementsRegions[r]->elements->datatarray().end(), elements->offset->datatarray().distribution()[t]);
			auto end = std::lower_bound(elementsRegions[r]->elements->datatarray().begin(), elementsRegions[r]->elements->datatarray().end(), elements->offset->datatarray().distribution()[t + 1]);
			for (auto i = begin; i != end; ++i) {
				regionElementMask[*i * eregionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto offset = elements->offset->datatarray().data();
		auto epointer = elements->epointers->datatarray().data();
		auto enodes = elements->nodes->cbegin(t);
		auto eneighbors = elements->faceNeighbors->cbegin(t);
		auto output = elements->inputOffset->cbegin(t);
		auto nIDs = nodes->IDs->datatarray().data();

		std::vector<std::vector<esint> > tsElements(targets.size(), std::vector<esint>({ 0 }));

		std::vector<esint>  telemsOffset;
		std::vector<Element*> telemsEpointer;
		std::vector<esint>  telemsNodesDistribution;
		std::vector<esint>  telemsNodesData;
		std::vector<esint>  telemsNeighborsDistribution;
		std::vector<esint>  telemsNeighborsData;
		std::vector<esint>  telemsOutputDistribution;
		std::vector<esint>  telemsOutputData;
		std::vector<esint>  telemsRegions;
		if (t == 0) {
			telemsNodesDistribution.push_back(0);
			telemsNeighborsDistribution.push_back(0);
			telemsOutputDistribution.push_back(0);
		}

		// estimation
		telemsOffset.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsEpointer.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsNodesDistribution.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsNeighborsDistribution.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsOutputDistribution.reserve(1.5 * elements->offset->datatarray().size() / threads);
		telemsRegions.reserve(1.5 * elements->offset->datatarray().size() / threads);

		size_t target;
		for (size_t e = elements->offset->datatarray().distribution()[t]; e < elements->offset->datatarray().distribution()[t + 1]; ++e, ++enodes, ++eneighbors, ++output) {
			PolyElement poly(Element::decode(epointer[e]->code, enodes->size()), enodes->begin());
			if (partition[e] == info::mpi::rank) {
				telemsOffset.push_back(offset[e]);
				telemsEpointer.push_back(epointer[e]);
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					if (poly.isNode(n - enodes->begin())) {
						telemsNodesData.push_back(nIDs[*n]);
					} else {
						telemsNodesData.push_back(*n);
					}
				}
				telemsNodesDistribution.push_back(telemsNodesData.size());
				telemsNeighborsData.insert(telemsNeighborsData.end(), eneighbors->begin(), eneighbors->end());
				telemsNeighborsDistribution.push_back(telemsNeighborsData.size());
				telemsOutputData.insert(telemsOutputData.end(), output->begin(), output->end());
				telemsOutputDistribution.push_back(telemsOutputData.size());
				telemsRegions.insert(telemsRegions.end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			} else {
				target = t2i(partition[e]);
				tsElements[target].insert(tsElements[target].end(), { offset[e], static_cast<int>(epointer[e]->code) });
				tsElements[target].push_back(enodes->size());
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					if (poly.isNode(n - enodes->begin())) {
						tsElements[target].push_back(nIDs[*n]);
					} else {
						tsElements[target].push_back(*n);
					}
				}
				tsElements[target].push_back(eneighbors->size());
				tsElements[target].insert(tsElements[target].end(), eneighbors->begin(), eneighbors->end());
				tsElements[target].push_back(output->size());
				tsElements[target].insert(tsElements[target].end(), output->begin(), output->end());
				tsElements[target].insert(tsElements[target].end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			}
		}

		elemsOffset[t].swap(telemsOffset);
		elemsEpointer[t].swap(telemsEpointer);
		elemsNodesDistribution[t].swap(telemsNodesDistribution);
		elemsNodesData[t].swap(telemsNodesData);
		elemsNeighborsDistribution[t].swap(telemsNeighborsDistribution);
		elemsNeighborsData[t].swap(telemsNeighborsData);
		elemsOutputDistribution[t].swap(telemsOutputDistribution);
		elemsOutputData[t].swap(telemsOutputData);
		elemsRegions[t].swap(telemsRegions);

		sElements[t].swap(tsElements);
	}

	profiler::synccheckpoint("serialize_elements");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE ELEMENTS");

	// Step 2: Serialize node data

	std::vector<esint> regionNodeMask(nodes->size * bregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->nodes) {
				maskOffset = r / (8 * sizeof(esint));
				esint bit = (esint)1 << (r % (8 * sizeof(esint)));
				auto begin = std::lower_bound(boundaryRegions[r]->nodes->datatarray().begin(), boundaryRegions[r]->nodes->datatarray().end(), nodes->IDs->datatarray().distribution()[t]);
				auto end = std::lower_bound(boundaryRegions[r]->nodes->datatarray().begin(), boundaryRegions[r]->nodes->datatarray().end(), nodes->IDs->datatarray().distribution()[t + 1]);
				for (auto i = begin; i != end; ++i) {
					regionNodeMask[*i * bregionsBitMaskSize + maskOffset] |= bit;
				}
			}
		}
	}

	esint eBegin = elements->offset->datatarray().size();
	Communication::exscan(eBegin);
	esint eEnd = eBegin + elements->offset->datatarray().size();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = nodes->IDs->datatarray();
		const auto &coordinates = nodes->coordinates->datatarray();
		auto elems = nodes->elements->cbegin(t);
		auto output = nodes->inputOffset->cbegin(t);

		std::vector<esint>  tnodesIDs;
		std::vector<Point>    tnodesCoordinates;
		std::vector<esint>  tnodesElemsDistribution;
		std::vector<esint>  tnodesElemsData;
		std::vector<esint>  tnodesOutputDistribution;
		std::vector<esint>  tnodesOutputData;
		std::vector<esint>  tnodesRegions;

		if (t == 0) {
			tnodesElemsDistribution.push_back(0);
			tnodesOutputDistribution.push_back(0);
		}

		std::vector<std::vector<esint> > tsNodes(targets.size(), std::vector<esint>({ 0 }));

		tnodesIDs.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesCoordinates.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesElemsDistribution.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesOutputDistribution.reserve(1.5 * nodes->IDs->datatarray().size() / threads);
		tnodesRegions.reserve(1.5 * nodes->IDs->datatarray().size() / threads);

		size_t target;
		std::vector<bool> last(targets.size() + 1); // targets + me
		for (size_t n = nodes->IDs->datatarray().distribution()[t]; n < nodes->IDs->datatarray().distribution()[t + 1]; ++n, ++elems, ++output) {
			std::fill(last.begin(), last.end(), false);
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				if (eBegin <= *e && *e < eEnd) {
					target = t2i(partition[*e - eBegin]);
					if (!last[target] && partition[*e - eBegin] != info::mpi::rank) {
						tsNodes[target].push_back(IDs[n]);
						tsNodes[target].insert(tsNodes[target].end(), reinterpret_cast<const esint*>(coordinates.data() + n), reinterpret_cast<const esint*>(coordinates.data() + n + 1));
						tsNodes[target].push_back(elems->size());
						tsNodes[target].insert(tsNodes[target].end(), elems->begin(), elems->end());
						tsNodes[target].push_back(output->size());
						tsNodes[target].insert(tsNodes[target].end(), output->begin(), output->end());
						tsNodes[target].insert(tsNodes[target].end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last[target] = true;
					}
					if (!last.back() && partition[*e - eBegin] == info::mpi::rank) {
						tnodesIDs.push_back(IDs[n]);
						tnodesCoordinates.push_back(coordinates[n]);
						tnodesElemsData.insert(tnodesElemsData.end(), elems->begin(), elems->end());
						tnodesElemsDistribution.push_back(tnodesElemsData.size());
						tnodesOutputData.insert(tnodesOutputData.end(), output->begin(), output->end());
						tnodesOutputDistribution.push_back(tnodesOutputData.size());
						tnodesRegions.insert(tnodesRegions.end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last.back() = true;
					}
				}
			}
		}

		nodesIDs[t].swap(tnodesIDs);
		nodesCoordinates[t].swap(tnodesCoordinates);
		nodesElemsDistribution[t].swap(tnodesElemsDistribution);
		nodesElemsData[t].swap(tnodesElemsData);
		nodesOutputDistribution[t].swap(tnodesOutputDistribution);
		nodesOutputData[t].swap(tnodesOutputData);
		nodesRegions[t].swap(tnodesRegions);

		sNodes[t].swap(tsNodes);
	}

	profiler::synccheckpoint("serialize_nodes");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE NODES");

	// Step 2.1: Serialize boundary regions data

	std::vector<esint> emembership;

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->originalDimension) {
			emembership.clear();
			emembership.resize(boundaryRegions[r]->epointers->datatarray().size());
			const std::vector<size_t> &distribution = boundaryRegions[r]->epointers->datatarray().distribution();

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto enodes = boundaryRegions[r]->elements->cbegin() + distribution[t];
				auto ep = boundaryRegions[r]->epointers->datatarray().begin(t);
				std::vector<esint> nlinks;
				int counter;
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes, ++ep) {
					PolyElement poly(Element::decode((*ep)->code, enodes->size()), enodes->begin());
					nlinks.clear();
					for (auto n = enodes->begin(); n != enodes->end(); ++n) {
						if (poly.isNode(n - enodes->begin())) {
							auto links = nodes->elements->cbegin() + *n;
							nlinks.insert(nlinks.end(), links->begin(), links->end());
						}
					}
					std::sort(nlinks.begin(), nlinks.end());
					counter = 1;
					for (size_t i = 1; i < nlinks.size(); ++i) {
						if (nlinks[i - 1] == nlinks[i]) {
							++counter;
							if (counter == poly.size && eBegin <= nlinks[i]) {
								emembership[e] = nlinks[i] - eBegin;
								break;
							}
						} else {
							counter = 1;
						}
					}
				}
			}

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				std::vector<esint>  tboundaryEDistribution;
				std::vector<esint>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;
				if (t == 0) {
					tboundaryEDistribution.push_back(0);
				}

				std::vector<size_t> tsize(targets.size());
				for (size_t i = 0; i < targets.size(); i++) {
					tsize[i] = sBoundary[t][i].size();
				}
				esint target;
				auto enodes = boundaryRegions[r]->elements->cbegin() + distribution[t];
				const auto &IDs = nodes->IDs->datatarray();
				const auto &epointer = boundaryRegions[r]->epointers->datatarray();
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					PolyElement poly(Element::decode(epointer[e]->code, enodes->size()), enodes->begin());
					if (partition[emembership[e]] == info::mpi::rank) {
						tboundaryEPointers.push_back(epointer[e]);
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							if (poly.isNode(n - enodes->begin())) {
								tboundaryEData.push_back(IDs[*n]);
							} else {
								tboundaryEData.push_back(*n);
							}
						}
						tboundaryEDistribution.push_back(tboundaryEData.size());
					} else {
						target = t2i(partition[emembership[e]]);
						sBoundary[t][target].push_back(static_cast<int>(epointer[e]->code));
						sBoundary[t][target].push_back(enodes->size());
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							if (poly.isNode(n - enodes->begin())) {
								sBoundary[t][target].push_back(IDs[*n]);
							} else {
								sBoundary[t][target].push_back(*n);
							}
						}
					}
				}
				for (size_t i = 0; i < targets.size(); i++) {
					sBoundary[t][i][r] = sBoundary[t][i].size() - tsize[i];
				}

				boundaryEDistribution[r][t].swap(tboundaryEDistribution);
				boundaryEData[r][t].swap(tboundaryEData);
				boundaryEPointers[r][t].swap(tboundaryEPointers);
			}
		}
	}

	profiler::synccheckpoint("serialize_boundary");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE BOUNDARIES");

	// Step 3: Send data to target processes

	#pragma omp parallel for
	for (size_t target = 0; target < targets.size(); ++target) {
		sElements[0][target].front() = sElements[0][target].size();
		sNodes[0][target].front() = sNodes[0][target].size();
		for (size_t t = 1; t < threads; t++) {
			sElements[t][target].front() = sElements[t][target].size();
			sElements[0][target].insert(sElements[0][target].end(), sElements[t][target].begin(), sElements[t][target].end());
			sNodes[t][target].front() = sNodes[t][target].size();
			sNodes[0][target].insert(sNodes[0][target].end(), sNodes[t][target].begin(), sNodes[t][target].end());
			sBoundary[0][target].insert(sBoundary[0][target].end(), sBoundary[t][target].begin(), sBoundary[t][target].end());
		}
	}

	if (!Communication::sendVariousTargets(sElements[0], rElements, targets)) {
		eslog::internalFailure("exchange elements data.\n");
	}

	if (!Communication::sendVariousTargets(sNodes[0], rNodes, targets)) {
		eslog::internalFailure("exchange nodes data.\n");
	}

	if (!Communication::sendVariousTargets(sBoundary[0], rBoundary, targets)) {
		eslog::internalFailure("exchange boundary data.\n");
	}

	profiler::synccheckpoint("exchange");
	eslog::checkpointln("EXCHANGE EL: EXCHANGE");

	// Step 4: Deserialize element data
	for (size_t i = 0; i < rElements.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rElements[i].size()) {
			rdistribution.push_back(p += rElements[i][p]);
		}

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<esint>  telemsOffset;
			std::vector<Element*> telemsEpointer;
			std::vector<esint>  telemsNodesDistribution;
			std::vector<esint>  telemsNodesData;
			std::vector<esint>  telemsNeighborsDistribution;
			std::vector<esint>  telemsNeighborsData;
			std::vector<esint>  telemsOutputDistribution;
			std::vector<esint>  telemsOutputData;
			std::vector<esint>  telemsRegions;

			telemsOffset.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsEpointer.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsNodesDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsNeighborsDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsOutputDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsRegions.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0, neighOffset = 0, outputOffset = 0;

			if (elemsNodesDistribution[t].size()) {
				distOffset = elemsNodesDistribution[t].back();
			}
			if (elemsNeighborsDistribution[t].size()) {
				neighOffset = elemsNeighborsDistribution[t].back();
			}
			if (elemsOutputDistribution[t].size()) {
				outputOffset = elemsOutputDistribution[t].back();
			}
			if (t == 0 && elemsNodesDistribution[t].size() == 0) {
				telemsNodesDistribution.push_back(0);
			}
			if (t == 0 && elemsNeighborsDistribution[t].size() == 0) {
				telemsNeighborsDistribution.push_back(0);
			}
			if (t == 0 && elemsOutputDistribution[t].size() == 0) {
				telemsOutputDistribution.push_back(0);
			}

			for (size_t e = rdistribution[t] + 1; e < rdistribution[t + 1]; ) {
				telemsOffset.push_back(rElements[i][e++]);
				telemsEpointer.push_back(&Mesh::edata[rElements[i][e++]]);
				telemsNodesData.insert(telemsNodesData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNodesDistribution.push_back(telemsNodesData.size() + distOffset);
				e += rElements[i][e++]; // nodes + nodes size
				telemsNeighborsData.insert(telemsNeighborsData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNeighborsDistribution.push_back(telemsNeighborsData.size() + neighOffset);
				e += rElements[i][e++]; // neighbors + neighbors size
				telemsOutputData.insert(telemsOutputData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsOutputDistribution.push_back(telemsOutputData.size() + outputOffset);
				e += rElements[i][e++]; // output + output size
				telemsRegions.insert(telemsRegions.end(), rElements[i].begin() + e, rElements[i].begin() + e + eregionsBitMaskSize);
				e += eregionsBitMaskSize;
			}

			elemsOffset[t].insert(elemsOffset[t].end(), telemsOffset.begin(), telemsOffset.end());
			elemsEpointer[t].insert(elemsEpointer[t].end(), telemsEpointer.begin(), telemsEpointer.end());
			elemsNodesDistribution[t].insert(elemsNodesDistribution[t].end(), telemsNodesDistribution.begin(), telemsNodesDistribution.end());
			elemsNodesData[t].insert(elemsNodesData[t].end(), telemsNodesData.begin(), telemsNodesData.end());
			elemsNeighborsDistribution[t].insert(elemsNeighborsDistribution[t].end(), telemsNeighborsDistribution.begin(), telemsNeighborsDistribution.end());
			elemsNeighborsData[t].insert(elemsNeighborsData[t].end(), telemsNeighborsData.begin(), telemsNeighborsData.end());
			elemsOutputDistribution[t].insert(elemsOutputDistribution[t].end(), telemsOutputDistribution.begin(), telemsOutputDistribution.end());
			elemsOutputData[t].insert(elemsOutputData[t].end(), telemsOutputData.begin(), telemsOutputData.end());
			elemsRegions[t].insert(elemsRegions[t].end(), telemsRegions.begin(), telemsRegions.end());
		}
	}

	profiler::synccheckpoint("deserialize_elements");
	eslog::checkpointln("EXCHANGE EL: DESERIALIZE ELEMENTS");

	// Step 4: Deserialize node data
	std::vector<esint> nodeset;
	for (size_t t = 0; t < threads; t++) {
		nodeset.insert(nodeset.end(), nodesIDs[t].begin(), nodesIDs[t].end());
	}
	std::sort(nodeset.begin(), nodeset.end());
	for (size_t i = 0; i < rNodes.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rNodes[i].size()) {
			rdistribution.push_back(p += rNodes[i][p]);
		}
		std::vector<std::vector<esint> > tnodeset(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point point;

			std::vector<esint>  tnodesIDs;
			std::vector<Point>    tnodesCoordinates;
			std::vector<esint>  tnodesElemsDistribution;
			std::vector<esint>  tnodesElemsData;
			std::vector<esint>  tnodesOutputDistribution;
			std::vector<esint>  tnodesOutputData;
			std::vector<esint>  tnodesRegions;
			std::vector<esint>  tnodeSet;
			std::vector<esint>  tnpermutation;

			tnodesIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesCoordinates.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesElemsDistribution.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesOutputDistribution.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesRegions.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodeSet.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0, outputOffset = 0;

			if (nodesElemsDistribution[t].size()) {
				distOffset = nodesElemsDistribution[t].back();
			}
			if (nodesOutputDistribution[t].size()) {
				outputOffset = nodesOutputDistribution[t].back();
			}
			if (t == 0 && nodesElemsDistribution[t].size() == 0) {
				tnodesElemsDistribution.push_back(0);
			}
			if (t == 0 && nodesOutputDistribution[t].size() == 0) {
				tnodesOutputDistribution.push_back(0);
			}

			for (size_t n = rdistribution[t] + 1; n < rdistribution[t + 1]; ) {
				tnpermutation.push_back(n);
				n += 1 + sizeof(Point) / sizeof(esint); // id, Point
				n += 1 + rNodes[i][n]; // linksize, links
				n += 1 + rNodes[i][n]; // outputsize, outputs
				n += bregionsBitMaskSize; // region mask
			}
			std::sort(tnpermutation.begin(), tnpermutation.end(), [&] (esint n1, esint n2) {
				return rNodes[i][n1] < rNodes[i][n2];
			});
			auto nodesetit = nodeset.begin();
			if (rdistribution[t] + 1 < rdistribution[t + 1]) {
				nodesetit = std::lower_bound(nodeset.begin(), nodeset.end(), rNodes[i][tnpermutation.front()]);
			}
			size_t index;
			for (size_t n = 0; n < tnpermutation.size(); n++) {
				index = tnpermutation[n];
				while (nodesetit != nodeset.end() && *nodesetit < rNodes[i][index]) ++nodesetit;
				if (nodesetit == nodeset.end() || *nodesetit != rNodes[i][index]) {
					tnodesIDs.push_back(rNodes[i][index]);
					tnodeSet.push_back(rNodes[i][index]);
					index += 1; //ID
					memcpy(reinterpret_cast<void*>(&point), rNodes[i].data() + index, sizeof(Point));
					tnodesCoordinates.push_back(point);
					index += sizeof(Point) / sizeof(esint); // points

					tnodesElemsData.insert(tnodesElemsData.end(), rNodes[i].begin() + index + 1, rNodes[i].begin() + index + 1 + rNodes[i][index]);
					tnodesElemsDistribution.push_back(tnodesElemsData.size() + distOffset);
					index += rNodes[i][index] + 1; // linksize + links

					tnodesOutputData.insert(tnodesOutputData.end(), rNodes[i].begin() + index + 1, rNodes[i].begin() + index + 1 + rNodes[i][index]);
					tnodesOutputDistribution.push_back(tnodesOutputData.size() + outputOffset);
					index += rNodes[i][index] + 1; // outputsize + output

					tnodesRegions.insert(tnodesRegions.end(), rNodes[i].begin() + index, rNodes[i].begin() + index + bregionsBitMaskSize);
					index += bregionsBitMaskSize; // region mask
				}
			}

			nodesIDs[t].insert(nodesIDs[t].end(), tnodesIDs.begin(), tnodesIDs.end());
			nodesCoordinates[t].insert(nodesCoordinates[t].end(), tnodesCoordinates.begin(), tnodesCoordinates.end());
			nodesElemsDistribution[t].insert(nodesElemsDistribution[t].end(), tnodesElemsDistribution.begin(), tnodesElemsDistribution.end());
			nodesElemsData[t].insert(nodesElemsData[t].end(), tnodesElemsData.begin(), tnodesElemsData.end());
			nodesOutputDistribution[t].insert(nodesOutputDistribution[t].end(), tnodesOutputDistribution.begin(), tnodesOutputDistribution.end());
			nodesOutputData[t].insert(nodesOutputData[t].end(), tnodesOutputData.begin(), tnodesOutputData.end());
			nodesRegions[t].insert(nodesRegions[t].end(), tnodesRegions.begin(), tnodesRegions.end());
			tnodeset[t].swap(tnodeSet);
		}

		size_t nsize = nodeset.size();
		for (size_t t = 0; t < threads; t++) {
			nodeset.insert(nodeset.end(), tnodeset[t].begin(), tnodeset[t].end());
		}
		std::inplace_merge(nodeset.begin(), nodeset.begin() + nsize, nodeset.end());
	}

	utils::threadDistributionToFullDistribution(elemsNodesDistribution);
	utils::threadDistributionToFullDistribution(elemsNeighborsDistribution);
	utils::threadDistributionToFullDistribution(elemsOutputDistribution);
	utils::threadDistributionToFullDistribution(nodesElemsDistribution);
	utils::threadDistributionToFullDistribution(nodesOutputDistribution);

	profiler::synccheckpoint("deserialize_nodes");
	eslog::checkpointln("EXCHANGE EL: DESERIALIZE NODES");

	// Step 4: Deserialize boundary data
	for (size_t n = 0; n < rBoundary.size(); ++n) {
		std::vector<std::vector<esint> > toffset(boundaryRegions.size()), tsize(boundaryRegions.size());
		esint offset = 0, p = 0;
		for (size_t t = 0; t < threads; t++) {
			for (size_t r = 0; r < boundaryRegions.size(); r++) {
				toffset[r].push_back(offset + boundaryRegions.size());
				tsize[r].push_back(rBoundary[n][p + r]);
				offset += rBoundary[n][p + r];
			}
			offset += boundaryRegions.size();
			p = offset;
		}
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {

				std::vector<esint>  tboundaryEDistribution;
				std::vector<esint>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;

				esint distOffset = 0;

				if (boundaryEDistribution[r][t].size()) {
					distOffset = boundaryEDistribution[r][t].back();
				}
				if (t == 0 && boundaryEDistribution[r][t].size() == 0) {
					tboundaryEDistribution.push_back(0);
				}

				for (esint i = toffset[r][t]; i < toffset[r][t] + tsize[r][t];) {
					tboundaryEPointers.push_back(&Mesh::edata[rBoundary[n][i++]]);
					tboundaryEData.insert(tboundaryEData.end(), rBoundary[n].begin() + i + 1, rBoundary[n].begin() + i + 1 + rBoundary[n][i]);
					tboundaryEDistribution.push_back(tboundaryEData.size() + distOffset);
					i += rBoundary[n][i] + 1;
				}

				boundaryEPointers[r][t].insert(boundaryEPointers[r][t].end(), tboundaryEPointers.begin(), tboundaryEPointers.end());
				boundaryEDistribution[r][t].insert(boundaryEDistribution[r][t].end(), tboundaryEDistribution.begin(), tboundaryEDistribution.end());
				boundaryEData[r][t].insert(boundaryEData[r][t].end(), tboundaryEData.begin(), tboundaryEData.end());
			}
		}
	}

	profiler::synccheckpoint("deserialize_boundary");
	eslog::checkpointln("EXCHANGE EL: DESERIALIZE BOUDARIES");

	// elements are redistributed later while decomposition -> distribution is not changed now
	std::vector<size_t> elemDistribution(threads + 1);
	for (size_t t = 1; t <= threads; t++) {
		elemDistribution[t] = elemDistribution[t - 1] + elemsOffset[t - 1].size();
	}

	newElements->offset = new serializededata<esint, esint>(1, elemsOffset);
	newElements->epointers = new serializededata<esint, Element*>(1, elemsEpointer);
	newElements->nodes = new serializededata<esint, esint>(elemsNodesDistribution, elemsNodesData); // global IDs

	newElements->regions = new serializededata<esint, esint>(eregionsBitMaskSize, elemsRegions);

	newElements->faceNeighbors = new serializededata<esint, esint>(elemsNeighborsDistribution, elemsNeighborsData);
	newElements->inputOffset = new serializededata<esint, esint>(elemsOutputDistribution, elemsOutputData);

	newElements->distribution.process.size = newElements->offset->structures();
	newElements->distribution.threads = newElements->offset->datatarray().distribution();

	// Step 5: Balance node data to threads
	std::vector<size_t> nodeDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		nodeDistribution[t] = nodeDistribution[t - 1] + nodesIDs[t - 1].size();
	}

	serializededata<esint, esint>::balance(1, nodesIDs);
	serializededata<esint, Point>::balance(1, nodesCoordinates);
	serializededata<esint, esint>::balance(nodesElemsDistribution, nodesElemsData);

	newNodes->IDs = new serializededata<esint, esint>(1, nodesIDs);
	newNodes->coordinates = new serializededata<esint, Point>(1, nodesCoordinates);
	newNodes->elements = new serializededata<esint, esint>(nodesElemsDistribution, nodesElemsData);
	newNodes->inputOffset = new serializededata<esint, esint>(nodesOutputDistribution, nodesOutputData);
	newNodes->size = newNodes->IDs->datatarray().size();
	newNodes->distribution = newNodes->IDs->datatarray().distribution();

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->originalDimension) {
			delete boundaryRegions[r]->elements;
			delete boundaryRegions[r]->epointers;

			utils::threadDistributionToFullDistribution(boundaryEDistribution[r]);
			boundaryRegions[r]->elements = new serializededata<esint, esint>(boundaryEDistribution[r], boundaryEData[r]);
			boundaryRegions[r]->epointers = new serializededata<esint, Element*>(1, boundaryEPointers[r]);
		}
	}

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		esint maskOffset = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));
		delete elementsRegions[r]->elements;
		std::vector<std::vector<esint> > regionelems(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = 0; i < elemsRegions[t].size(); i += eregionsBitMaskSize) {
				if (elemsRegions[t][i + maskOffset] & bit) {
					regionelems[t].push_back(elemDistribution[t] + i / eregionsBitMaskSize);
				}
			}
		}

		serializededata<esint, esint>::balance(1, regionelems);
		elementsRegions[r]->elements = new serializededata<esint, esint>(1, regionelems);
	}

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->nodes) {
			esint maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));
			delete boundaryRegions[r]->nodes;
			std::vector<std::vector<esint> > regionnodes(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = 0; i < nodesRegions[t].size(); i += bregionsBitMaskSize) {
					if (nodesRegions[t][i + maskOffset] & bit) {
						regionnodes[t].push_back(nodeDistribution[t] + i / bregionsBitMaskSize);
					}
				}
			}

			serializededata<esint, esint>::balance(1, regionnodes);
			boundaryRegions[r]->nodes = new serializededata<esint, esint>(1, regionnodes);
		}
	}

	std::vector<size_t> eIDsOLD = Communication::getDistribution(elements->epointers->datatarray().size());
	std::vector<size_t> eIDsNEW = Communication::getDistribution(newElements->epointers->datatarray().size());

	for (size_t t = 1; t < threads; ++t) {
		elemsOffset[0].insert(elemsOffset[0].end(), elemsOffset[t].begin(), elemsOffset[t].end());
	}

	std::vector<esint> epermutation(elemsOffset[0].size());
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return elemsOffset[0][i] < elemsOffset[0][j]; });

	std::vector<esint> sortedElements(epermutation.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		for (size_t e = elemDistribution[t]; e < elemDistribution[t + 1]; e++) {
			sortedElements[e] = elemsOffset[0][epermutation[e]];
		}
		utils::sortAndRemoveDuplicates(nodesElemsData[t]);
	}
	utils::inplaceMerge(nodesElemsData);
	utils::removeDuplicates(nodesElemsData[0]);

	std::vector<std::vector<esint> > requestedIDs, receivedTargets, IDrequests;
	std::vector<std::vector<int> > IDtargets(threads);
	std::vector<int> sources;

	std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, info::mpi::size);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		std::vector<int> ttargets;

		auto eIDbegin = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[rdistribution[t]]);
		auto eIDend   = std::lower_bound(eIDbegin, nodesElemsData[0].end(), eIDsOLD[rdistribution[t + 1]]);

		for (size_t r = rdistribution[t]; r < rdistribution[t + 1]; r++) {
			eIDbegin = std::lower_bound(eIDbegin, eIDend, eIDsOLD[r]);
			auto end = std::lower_bound(eIDbegin, eIDend, eIDsOLD[r + 1]);
			if (eIDbegin != end) {
				ttargets.push_back(r);
			}
			eIDbegin = end;
		}

		IDtargets[t].swap(ttargets);
	}

	for (size_t t = 1; t < threads; ++t) {
		IDtargets[0].insert(IDtargets[0].end(), IDtargets[t].begin(), IDtargets[t].end());
	}

	requestedIDs.resize(IDtargets[0].size());

	#pragma omp parallel for
	for (size_t t = 0; t < IDtargets[0].size(); t++) {
		auto mybegin = std::lower_bound(sortedElements.begin(), sortedElements.end(), eIDsOLD[IDtargets[0][t]]);
		auto myend   = std::lower_bound(sortedElements.begin(), sortedElements.end(), eIDsOLD[IDtargets[0][t] + 1]);
		auto nbegin = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[IDtargets[0][t]]);
		auto nend   = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[IDtargets[0][t] + 1]);
		requestedIDs[t].resize(nend - nbegin);
		requestedIDs[t].resize(std::set_difference(nbegin, nend, mybegin, myend, requestedIDs[t].begin()) - requestedIDs[t].begin());
	}

	profiler::synccheckpoint("postprocessing");
	eslog::checkpointln("EXCHANGE EL: POST-PROCESSING");

	if (!Communication::sendVariousTargets(requestedIDs, IDrequests, IDtargets[0], sources)) {
		eslog::internalFailure("exchange ID requests.\n");
	}

	for (size_t r = 0; r < IDrequests.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, IDrequests[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			if (rdistribution[t] != rdistribution[t + 1]) {
				for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
					IDrequests[r][e] = partition[IDrequests[r][e] - eIDsOLD[info::mpi::rank]];
				}
			}
		}
	}

	if (!Communication::sendVariousTargets(IDrequests, receivedTargets, sources)) {
		eslog::internalFailure("return ID targets.\n");
	}

	profiler::synccheckpoint("idrequests");
	eslog::checkpointln("EXCHANGE EL: IDS REQUESTS");

	IDtargets.clear();
	IDtargets.resize(threads);

	std::vector<std::vector<bool> > newtargets(threads, std::vector<bool>(info::mpi::size, false));
	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				newtargets[t][receivedTargets[r][e]] = true;
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		std::vector<int> ttargets;
		for (int r = 0; r < info::mpi::size; r++) {
			if (r != info::mpi::rank && newtargets[t][r]) {
				ttargets.push_back(r);
			}
		}
		IDtargets[t].swap(ttargets);
	}
	for (size_t t = 1; t < threads; ++t) {
		IDtargets[0].insert(IDtargets[0].end(), IDtargets[t].begin(), IDtargets[t].end());
	}
	utils::sortAndRemoveDuplicates(IDtargets[0]);

	std::vector<std::vector<std::vector<esint> > > newIDrequests(threads, std::vector<std::vector<esint> >(IDtargets[0].size()));
	std::vector<std::vector<esint> > newIDs;

	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			size_t tindex;
			std::vector<std::vector<esint> > tnewIDrequests(IDtargets[0].size());

			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				tindex = std::lower_bound(IDtargets[0].begin(), IDtargets[0].end(), receivedTargets[r][e]) - IDtargets[0].begin();
				tnewIDrequests[tindex].push_back(requestedIDs[r][e]);
			}

			for (size_t n = 0; n < IDtargets[0].size(); n++) {
				newIDrequests[t][n].insert(newIDrequests[t][n].end(), tnewIDrequests[n].begin(), tnewIDrequests[n].end());
			}
		}

		for (size_t t = 1; t < threads; ++t) {
			for (size_t n = 0; n < IDtargets[0].size(); n++) {
				newIDrequests[0][n].insert(newIDrequests[0][n].end(), newIDrequests[t][n].begin(), newIDrequests[t][n].end());
				newIDrequests[t][n].clear();
			}
		}
	}

	profiler::synccheckpoint("process_requests");
	eslog::checkpointln("EXCHANGE EL: PROC REQUESTS");

	if (!Communication::sendVariousTargets(newIDrequests[0], IDrequests, IDtargets[0], sources)) {
		eslog::internalFailure("request new ID.\n");
	}

	for (size_t r = 0; r < IDrequests.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, IDrequests[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			if (rdistribution[t] != rdistribution[t + 1]) {
				auto eit = std::lower_bound(sortedElements.begin(), sortedElements.end(), IDrequests[r][rdistribution[t]]);
				for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
					while (eit != sortedElements.end() && *eit < IDrequests[r][e]) ++eit;
					IDrequests[r][e] = epermutation[eit - sortedElements.begin()] + eIDsNEW[info::mpi::rank];
				}
			}
		}
	}

	if (!Communication::sendVariousTargets(IDrequests, newIDs, sources)) {
		eslog::internalFailure("return new ID.\n");
	}

	profiler::synccheckpoint("exchange_ids");
	eslog::checkpointln("EXCHANGE EL: NEW IDS");

	std::vector<size_t> offsets = { 0 };
	for (size_t r = 0; r < newIDs.size(); r++) {
		offsets.push_back(offsets.back() + newIDs[r].size());
	}
	std::vector<std::pair<esint, esint> > IDMap(offsets.back());

	for (size_t r = 0; r < newIDs.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, newIDrequests[0][r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				IDMap[offsets[r] + e].first = newIDrequests[0][r][e];
				IDMap[offsets[r] + e].second = newIDs[r][e];
			}
		}
	}

	utils::sortWithInplaceMerge(IDMap, offsets);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto elem = newNodes->elements->begin(t); elem != newNodes->elements->end(t); ++elem) {
			for (auto e = elem->begin(); e != elem->end(); ++e) {
				auto mapit = std::lower_bound(IDMap.begin(), IDMap.end(), std::make_pair(*e, (esint)0));
				if (mapit == IDMap.end() || mapit->first != *e) {
					*e = epermutation[std::lower_bound(sortedElements.begin(), sortedElements.end(), *e) - sortedElements.begin()] + eIDsNEW[info::mpi::rank];
				} else {
					*e = mapit->second;
				}
			}
			std::sort(elem->begin(), elem->end());
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto neighbors = newElements->faceNeighbors->begin(t); neighbors != newElements->faceNeighbors->end(t); ++neighbors) {
			for (auto n = neighbors->begin(); n != neighbors->end(); ++n) {
				if (*n != -1) {
					auto mapit = std::lower_bound(IDMap.begin(), IDMap.end(), std::make_pair(*n, (esint)0));
					if (mapit == IDMap.end() || mapit->first != *n) {
						*n = epermutation[std::lower_bound(sortedElements.begin(), sortedElements.end(), *n) - sortedElements.begin()] + eIDsNEW[info::mpi::rank];
					} else {
						*n = mapit->second;
					}
				}
			}
		}
	}


	std::vector<esint> permutation(newNodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return newNodes->IDs->datatarray()[i] < newNodes->IDs->datatarray()[j]; });

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ep = newElements->epointers->datatarray().begin(t);
		for (auto enodes = newElements->nodes->begin(t); enodes != newElements->nodes->end(t); ++enodes, ++ep) {
			PolyElement poly(Element::decode((*ep)->code, enodes->size()), enodes->begin());
			for (size_t n = 0; n < enodes->size(); ++n) {
				if (poly.isNode(n)) {
					enodes->at(n) = *std::lower_bound(permutation.begin(), permutation.end(), enodes->at(n), [&] (esint i, esint val) {
						return newNodes->IDs->datatarray()[i] < val;
					});
				}
			}
		}
	}

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->elements) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto ep = boundaryRegions[r]->epointers->datatarray().begin(t);
				for (auto enodes = boundaryRegions[r]->elements->begin(t); enodes != boundaryRegions[r]->elements->end(t); ++enodes, ++ep) {
					PolyElement poly(Element::decode((*ep)->code, enodes->size()), enodes->begin());
					for (size_t n = 0; n < enodes->size(); ++n) {
						if (poly.isNode(n)) {
							enodes->at(n) = *std::lower_bound(permutation.begin(), permutation.end(), enodes->at(n), [&] (esint i, esint val) {
								return newNodes->IDs->datatarray()[i] < val;
							});
						}
					}
				}
			}
		}
	}

	std::vector<std::vector<esint> > rankBoundaries(threads);
	std::vector<std::vector<int> > rankData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> trankBoundaries;
		std::vector<int> trankData;
		if (t == 0) {
			trankBoundaries.push_back(0);
		}

		int rank;
		for (auto elem = newNodes->elements->begin(t); elem != newNodes->elements->end(t); ++elem) {
			trankData.push_back(std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *elem->begin() + 1) - eIDsNEW.begin() - 1);
			for (auto e = elem->begin() + 1; e != elem->end(); ++e) {
				rank = std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *e + 1) - eIDsNEW.begin() - 1;
				if (rank != trankData.back()) {
					trankData.push_back(rank);
				}
			}
			trankBoundaries.push_back(trankData.size());
		}

		rankBoundaries[t].swap(trankBoundaries);
		rankData[t].swap(trankData);
	}

	utils::threadDistributionToFullDistribution(rankBoundaries);

	newNodes->ranks = new serializededata<esint, int>(rankBoundaries, rankData);

	std::iota(newElements->offset->datatarray().begin(), newElements->offset->datatarray().end(), eIDsNEW[info::mpi::rank]);
	newElements->distribution.process.offset = eIDsNEW[info::mpi::rank];
	newElements->distribution.process.totalSize = elements->distribution.process.totalSize;
	std::swap(elements, newElements);
	std::swap(nodes, newNodes);
	neighbors.clear();
	for (size_t t = 0; t < IDtargets[0].size(); t++) {
		if (IDtargets[0][t] != info::mpi::rank) {
			neighbors.push_back(IDtargets[0][t]);
		}
	}
	neighborsWithMe = neighbors;
	neighborsWithMe.push_back(info::mpi::rank);
	std::sort(neighborsWithMe.begin(), neighborsWithMe.end());

	delete newElements;
	delete newNodes;

	profiler::synccheckpoint("finish");
	profiler::syncend("exchange_elements");
	eslog::endln("EXCHANGE EL: FINISH");
	eslog::checkpointln("MESH: ELEMENTS EXCHANGED");
}

void sortNodes(NodeStore *nodes, ElementStore *elements, std::vector<BoundaryRegionStore*> &boundaryRegions)
{
	profiler::syncstart("sort_nodes");
	profiler::syncparam("size", nodes->size);
	std::vector<esint> permutation(nodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		int irank = *(nodes->ranks->begin() + i)->begin();
		int jrank = *(nodes->ranks->begin() + j)->begin();
		if (irank == jrank) {
			return nodes->IDs->datatarray()[i] < nodes->IDs->datatarray()[j];
		}
		return irank < jrank;
	});

	nodes->permute(permutation);
	profiler::synccheckpoint("permute");

	// nhalo is used by other routines (before arrange boundary elements)
	nodes->uniqInfo.nhalo = 0;
	auto ranks = nodes->ranks->begin();
	while (ranks != nodes->ranks->end() && *ranks->begin() < info::mpi::rank) {
		++ranks;
		++nodes->uniqInfo.nhalo;
	}

	std::vector<esint> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return permutation[i] < permutation[j]; });

	profiler::synccheckpoint("backpermute");

	auto localremap = [&] (serializededata<esint, Element*>* epointers, serializededata<esint, esint>* data) {
		if (data == NULL) {
			return;
		}
		#pragma omp parallel for
		for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
			auto ep = epointers->datatarray().begin(t);
			for (auto e = data->begin(t); e != data->end(t); ++e, ++ep) {
				PolyElement poly(Element::decode((*ep)->code, e->size()), e->begin());
				for (auto n = e->begin(); n != e->end(); ++n) {
					if (poly.isNode(n - e->begin())) {
						*n = backpermutation[*n];
					}
				}
			}
		}
	};

	localremap(elements->epointers, elements->nodes);

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->elements != NULL) {
			localremap(boundaryRegions[r]->epointers, boundaryRegions[r]->elements);
		}
		if (!StringCompare::caseInsensitiveEq(boundaryRegions[r]->name, "ALL_NODES")) {
			if (boundaryRegions[r]->nodes != NULL) {
			#pragma omp parallel for
			for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
				for (auto n = boundaryRegions[r]->nodes->datatarray().begin(t); n != boundaryRegions[r]->nodes->datatarray().end(t); ++n) {
					*n = backpermutation[*n];
				}
			}
				std::sort(boundaryRegions[r]->nodes->datatarray().begin(), boundaryRegions[r]->nodes->datatarray().end());
			}
		}
	}

	profiler::syncend("sort_nodes");
	eslog::checkpointln("MESH: NODES SORTED");
}

void computeElementDistribution(ElementStore *elements)
{
	profiler::syncstart("compute_element_distribution");
	elements->distribution.clear();
	auto enodes = elements->nodes->cbegin();
	for (auto epointer = elements->epointers->datatarray().begin(); epointer != elements->epointers->datatarray().end(); ++epointer, ++enodes) {
		++elements->distribution.code[(int)(*epointer)->code].size;
		if ((*epointer)->code == Element::CODE::POLYGON) {
			elements->distribution.polygonNodes.size += enodes->size() - 1;
		}
		if ((*epointer)->code == Element::CODE::POLYHEDRON) {
			elements->distribution.polyhedronFaces.size += enodes->front();
			elements->distribution.polyhedronNodes.size += enodes->size() - enodes->front() - 1;
		}
	}
	elements->distribution.process.size = elements->epointers->datatarray().size();

	std::vector<esint> sum, offset;
	for (size_t i = 0; i < elements->distribution.code.size(); ++i) {
		offset.push_back(elements->distribution.code[i].size);
	}
	offset.push_back(elements->distribution.polygonNodes.size);
	offset.push_back(elements->distribution.polyhedronFaces.size);
	offset.push_back(elements->distribution.polyhedronNodes.size);

	sum.resize(offset.size());
	Communication::exscan(sum, offset);
	for (size_t i = 0; i < elements->distribution.code.size(); ++i) {
		elements->distribution.code[i].offset = offset[i];
		elements->distribution.code[i].totalSize = sum[i];
	}
	elements->distribution.polygonNodes.offset    = offset[elements->distribution.code.size() + 0];
	elements->distribution.polyhedronNodes.offset = offset[elements->distribution.code.size() + 1];
	elements->distribution.polyhedronFaces.offset = offset[elements->distribution.code.size() + 2];
	elements->distribution.polygonNodes.totalSize    = sum[elements->distribution.code.size() + 0];
	elements->distribution.polyhedronNodes.totalSize = sum[elements->distribution.code.size() + 1];
	elements->distribution.polyhedronFaces.totalSize = sum[elements->distribution.code.size() + 2];

	elements->distribution.process.offset = elements->distribution.process.size;
	elements->distribution.process.totalSize = Communication::exscan(elements->distribution.process.offset);
	elements->distribution.process.next = elements->distribution.process.offset + elements->distribution.process.size;
	profiler::syncend("compute_element_distribution");

	eslog::checkpointln("MESH: ELEMENT DISTRIBUTION COMPUTED");
}

void computeERegionsDistribution(ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions)
{
	profiler::syncstart("regions element distribution");
	std::vector<esint> sum, offset;
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		elementsRegions[r]->distribution.clear();
		auto enodes = elements->nodes->cbegin();
		for (auto e = elementsRegions[r]->elements->datatarray().begin(); e != elementsRegions[r]->elements->datatarray().end(); ++e, ++enodes) {
			++elementsRegions[r]->distribution.code[(int)elements->epointers->datatarray()[*e]->code].size;
			if (elements->epointers->datatarray()[*e]->code == Element::CODE::POLYGON) {
				elementsRegions[r]->distribution.polygonNodes.size += enodes->size() - 1;
			}
			if (elements->epointers->datatarray()[*e]->code == Element::CODE::POLYHEDRON) {
				elementsRegions[r]->distribution.polyhedronFaces.size += enodes->front();
				elementsRegions[r]->distribution.polyhedronNodes.size += enodes->size() - enodes->front() - 1;
			}
		}

		for (size_t i = 0; i < elements->distribution.code.size(); ++i) {
			elementsRegions[r]->distribution.process.size += elementsRegions[r]->distribution.code[i].size;
			elementsRegions[r]->distribution.process.next += elementsRegions[r]->distribution.code[i].size;
			offset.push_back(elementsRegions[r]->distribution.code[i].size);
			if (Mesh::element(i).code == Element::CODE::POLYGON) {
				offset.push_back(elementsRegions[r]->distribution.polygonNodes.size);
			}
			if (Mesh::element(i).code == Element::CODE::POLYHEDRON) {
				offset.push_back(elementsRegions[r]->distribution.polyhedronFaces.size);
				offset.push_back(elementsRegions[r]->distribution.polyhedronNodes.size);
			}
		}
		elementsRegions[r]->distribution.process.offset = elementsRegions[r]->distribution.process.size;
		offset.push_back(elementsRegions[r]->distribution.process.offset);
	}

	sum.resize(offset.size());
	Communication::exscan(sum, offset);

	for (size_t r = 0, j = 0; r < elementsRegions.size(); r++) {
		for (size_t i = 0; i < elements->distribution.code.size(); i++) {
			elementsRegions[r]->distribution.code[i].offset = offset[j];
			elementsRegions[r]->distribution.code[i].totalSize = sum[j++];
			if (Mesh::element(i).code == Element::CODE::POLYGON) {
				elementsRegions[r]->distribution.polygonNodes.offset = offset[j];
				elementsRegions[r]->distribution.polygonNodes.totalSize = sum[j++];
			}
			if (Mesh::element(i).code == Element::CODE::POLYHEDRON) {
				elementsRegions[r]->distribution.polyhedronFaces.offset = offset[j];
				elementsRegions[r]->distribution.polyhedronFaces.totalSize = sum[j++];
				elementsRegions[r]->distribution.polyhedronNodes.offset = offset[j];
				elementsRegions[r]->distribution.polyhedronNodes.totalSize = sum[j++];
			}
		}
		elementsRegions[r]->distribution.process.offset = offset[j];
		elementsRegions[r]->distribution.process.next += offset[j];
		elementsRegions[r]->distribution.process.totalSize = sum[j++];
	}
	profiler::syncend("regions element distribution");

	elements->distribution.process = elementsRegions.front()->distribution.process;
	elements->distribution.code = elementsRegions.front()->distribution.code;
	elements->distribution.polygonNodes = elementsRegions.front()->distribution.polygonNodes;
	elements->distribution.polyhedronNodes = elementsRegions.front()->distribution.polyhedronNodes;
	elements->distribution.polyhedronFaces = elementsRegions.front()->distribution.polyhedronFaces;

	eslog::checkpointln("MESH: EREGIONS DISTRIBUTION COMPUTED");
}

} // namespace mesh
} // namespace espreso
