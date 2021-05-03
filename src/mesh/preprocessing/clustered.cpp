
#include "meshpreprocessing.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/contactinterfacestore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/bodystore.h"
#include "wrappers/mpi/communication.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace espreso {
namespace mesh {

void fillRegionMask(ElementStore *elements, const std::vector<ElementsRegionStore*> &elementsRegions)
{
	profiler::syncstart("fill_region_mask");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > eregions(threads);

	// regions are transfered via mask
	auto bitMastSize = [](size_t elements)
	{
		size_t size = 8 * sizeof(esint);
		return elements / size + (elements % size ? 1 : 0);
	};
	int regionsBitMaskSize = bitMastSize(elementsRegions.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		eregions[t].resize(regionsBitMaskSize * elements->epointers->datatarray().size(t));
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));

			const auto &relements = elementsRegions[r]->elements->datatarray();
			auto begin = std::lower_bound(relements.begin(), relements.end(), elements->epointers->datatarray().distribution()[t]);
			auto end = std::lower_bound(relements.begin(), relements.end(), elements->epointers->datatarray().distribution()[t + 1]);
			for (auto i = begin; i != end; ++i) {
				eregions[t][(*i - elements->epointers->datatarray().distribution()[t]) * regionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}
	elements->regions = new serializededata<esint, esint>(regionsBitMaskSize, eregions);

	profiler::syncend("fill_region_mask");
	eslog::checkpointln("MESH: REGION MASK FILLED");
}

ElementStore* exchangeHalo(ElementStore *elements, NodeStore *nodes, std::vector<int> &neighbors)
{
	profiler::syncstart("exchange_halo");
	// halo elements are all elements that have some shared node

	esint ebegin = elements->distribution.process.offset;
	esint eend = ebegin + elements->distribution.process.size;

	size_t threads = info::env::OMP_NUM_THREADS;
	std::vector<std::vector<esint> > sBuffer(neighbors.size()), rBuffer(neighbors.size());

	std::vector<std::vector<std::vector<esint> > > hElements(threads);

	// we have to got through all nodes because intervals are not computed yet
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<esint> > telements(neighbors.size());
		auto elinks = nodes->elements->cbegin(t);
		size_t i = 0;

		for (auto ranks = nodes->ranks->cbegin(t); ranks != nodes->ranks->cend(t); ++ranks, ++elinks) {
			auto begin = elinks->begin();
			auto end = elinks->begin();
			if (ranks->size() > 1) {
				i = 0;
				while (begin != elinks->end() && *begin < ebegin) ++begin;
				end = begin;
				while (end != elinks->end() && *end < eend) ++end;
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != info::mpi::rank) {
						while (neighbors[i] < *rank) ++i;
						telements[i].insert(telements[i].end(), begin, end);
					}
				}
			}
		}
		hElements[t].swap(telements);
	}

	int rsize = elements->regions->edataSize();

	std::vector<std::vector<size_t> > tdist(neighbors.size());
	for (size_t n = 0; n < neighbors.size(); ++n) {
		tdist[n] = { 0, hElements[0][n].size() };
	}
	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < neighbors.size(); ++n) {
			hElements[0][n].insert(hElements[0][n].end(), hElements[t][n].begin(), hElements[t][n].end());
			tdist[n].push_back(hElements[0][n].size());
		}
	}
	for (size_t n = 0; n < neighbors.size(); ++n) {
		utils::sortWithInplaceMerge(hElements[0][n], tdist[n]);
	}
	#pragma omp parallel for
	for (size_t n = 0; n < neighbors.size(); ++n) {
		utils::removeDuplicates(hElements[0][n]);
		tdist[n] = tarray<size_t>::distribute(threads, hElements[0][n].size());
		sBuffer[n].resize((4 + rsize) * hElements[0][n].size());
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = elements->IDs->datatarray();
		const auto &body = elements->body->datatarray();
		const auto &material = elements->material->datatarray();
		const auto &epointer = elements->epointers->datatarray();
		const auto &regions = elements->regions->datatarray();
		for (size_t n = 0; n < neighbors.size(); ++n) {
			for (size_t e = tdist[n][t]; e < tdist[n][t + 1]; e++) {
				sBuffer[n][(4 + rsize) * e + 0] = IDs[hElements[0][n][e] - ebegin];
				sBuffer[n][(4 + rsize) * e + 1] = body[hElements[0][n][e] - ebegin];
				sBuffer[n][(4 + rsize) * e + 2] = material[hElements[0][n][e] - ebegin];
				sBuffer[n][(4 + rsize) * e + 3] = (esint)epointer[hElements[0][n][e] - ebegin]->code;
				memcpy(sBuffer[n].data() + (4 + rsize) * e + 4, regions.data() + (hElements[0][n][e] - ebegin) * rsize, sizeof(esint) * rsize);
			}
		}
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighbors)) {
		eslog::internalFailure("exchange halo elements.\n");
	}
	profiler::synccheckpoint("exchange");

	std::vector<std::vector<esint> > hid(threads), hregions(threads);
	std::vector<std::vector<int> > hbody(threads), hmaterial(threads);

	std::vector<std::vector<Element*> > hcode(threads);

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, rBuffer[n].size() / (4 + rsize));
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; ++e) {
				hid[t].push_back(rBuffer[n][(4 + rsize) * e + 0]);
				hbody[t].push_back(rBuffer[n][(4 + rsize) * e + 1]);
				hmaterial[t].push_back(rBuffer[n][(4 + rsize) * e + 2]);
				hcode[t].push_back(&Mesh::edata[rBuffer[n][(4 + rsize) * e + 3]]);
				hregions[t].insert(hregions[t].end(), rBuffer[n].data() + (4 + rsize) * e + 4, rBuffer[n].data() + (4 + rsize) * e + 4 + rsize);
			}
		}
	}

	ElementStore *halo = new ElementStore();
	halo->IDs = new serializededata<esint, esint>(1, hid);
	halo->body = new serializededata<esint, int>(1, hbody);
	halo->material = new serializededata<esint, int>(1, hmaterial);
	halo->epointers = new serializededata<esint, Element*>(1, hcode);
	halo->regions = new serializededata<esint, esint>(rsize, hregions);

	halo->distribution.process.size = halo->IDs->datatarray().size();
	halo->distribution.threads = halo->IDs->datatarray().distribution();

	const auto &hIDs = halo->IDs->datatarray();
	std::vector<esint> permutation(hIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return hIDs[i] < hIDs[j]; });
	halo->permute(permutation);

	profiler::synccheckpoint("rbuffer");
	profiler::syncend("exchange_halo");
	eslog::checkpointln("MESH: HALO EXCHANGED");
	return halo;
}

void computeRegionsBoundaryElementsFromNodes(const NodeStore *nodes, const ElementStore *elements, const ElementStore *halo, const std::vector<ElementsRegionStore*> &elementsRegions, BoundaryRegionStore *bregion)
{
	profiler::syncstart("compute_boundary_elements_from_nodes");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<std::pair<esint, esint> > > elems(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = bregion->nodes->begin(t)->begin(); n != bregion->nodes->end(t)->begin(); ++n) {
			auto links = nodes->elements->begin() + *n;
			for (auto e = links->begin(); e != links->end(); ++e) {
				elems[t].push_back(std::make_pair(*e, *n));
			}
		}
	}

	std::vector<size_t> distribution = { 0, elems[0].size() };
	for (size_t t = 1; t < threads; t++) {
		elems[0].insert(elems[0].end(), elems[t].begin(), elems[t].end());
		distribution.push_back(elems[0].size());
	}

	utils::sortWithInplaceMerge(elems[0], distribution);

	auto begin = std::lower_bound(elems[0].begin(), elems[0].end(), elements->distribution.process.offset, [] (const std::pair<esint, esint> &p, esint e) { return p.first < e; });
	auto end = std::lower_bound(elems[0].begin(), elems[0].end(), elements->distribution.process.next, [] (const std::pair<esint, esint> &p, esint e) { return p.first < e; });

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);
	for (size_t t = 1; t < threads; t++) {
		while (
				begin + tdistribution[t] < end && begin <= begin + tdistribution[t] - 1 &&
				(begin + tdistribution[t] - 1)->first == (begin + tdistribution[t])->first) {

			++tdistribution[t];
		}
	}

	std::vector<std::vector<esint> > edist(threads), edata(threads), ecode(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	int rsize = elements->regions->edataSize();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> nodes, facenodes, lowerElements, lenodes;
		std::vector<esint> tdist, tdata, tcode;
		if (t == 0) {
			tdist.push_back(0);
		}

		int nface;
		esint element, neighbor, prev = 0;
		auto enodes = elements->nodes->cbegin();
		auto neighbors = elements->faceNeighbors->cbegin();
		const auto &regions = elements->regions->datatarray();

		for (size_t e = tdistribution[t]; e < tdistribution[t + 1]; e++) {
			nodes.push_back((begin + e)->second);
			if ((e + 1 == tdistribution[t + 1] || (begin + e + 1)->first != (begin + e)->first)) {

				element = (begin + e)->first - elements->distribution.process.offset;
				utils::sortAndRemoveDuplicates(nodes);

				enodes += element - prev;
				neighbors += element - prev;
				prev = element;

				const auto &fpointers = elements->epointers->datatarray()[element]->facepointers->datatarray();
				auto fnodes = elements->epointers->datatarray()[element]->faces->cbegin();
				nface = 0;
				for (auto f = fpointers.begin(); f != fpointers.end(); ++f, ++fnodes, ++nface) {

					auto addFace = [&] () {
						for (auto n = fnodes->begin(); n != fnodes->end(); ++n) {
							tdata.push_back(enodes->at(*n));
						}
						tdist.push_back(tdata.size());
						tcode.push_back((esint)(*f)->code);
					};

					if ((int)nodes.size() >= (*f)->nodes) {
						for (auto n = fnodes->begin(); n != fnodes->end(); ++n) {
							facenodes.push_back(enodes->at(*n));
						}
						std::sort(facenodes.begin(), facenodes.end());
						if (std::includes(nodes.begin(), nodes.end(), facenodes.begin(), facenodes.end())) {
							neighbor = neighbors->at(nface);
							if (neighbor == -1) {
								addFace();
							} else if (element + elements->distribution.process.offset < neighbor) {
								if (elements->distribution.process.isLocal(neighbor)) {
									neighbor -= elements->distribution.process.offset;
									if (memcmp(regions.data() + element * rsize, regions.data() + neighbor * rsize, sizeof(esint) * rsize) != 0) {
										addFace();
									}
								} else {
									neighbor = std::lower_bound(halo->IDs->datatarray().begin(), halo->IDs->datatarray().end(), neighbor) - halo->IDs->datatarray().begin();
									if (memcmp(regions.data() + element * rsize, halo->regions->datatarray().data() + neighbor * rsize, sizeof(esint) * rsize) != 0) {
										addFace();
									}
								}
							}
						}
						facenodes.clear();
					}
				}
				nodes.clear();
			}
		}

		edist[t].swap(tdist);
		edata[t].swap(tdata);
		ecode[t].swap(tcode);
	}

	utils::threadDistributionToFullDistribution(edist);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = 0; e < ecode[t].size(); e++) {
			epointers[t].push_back(&Mesh::edata[ecode[t][e]]);
		}
	}

	bregion->elements = new serializededata<esint, esint>(edist, edata);
	bregion->epointers = new serializededata<esint, Element*>(1, epointers);

	profiler::syncend("compute_boundary_elements_from_nodes");
	eslog::checkpoint("MESH: BOUNDARY FROM NODES COMPUTED");
	eslog::param("BOUNDARY", bregion->name.c_str());
	eslog::ln();
}

void computeRegionsBoundaryDistribution(std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces)
{
	profiler::syncstart("compute_regions_boudary_distribution");
	std::vector<BoundaryRegionStore*> allRegions;
	allRegions.insert(allRegions.end(), boundaryRegions.begin(), boundaryRegions.end());
	allRegions.insert(allRegions.end(), contactInterfaces.begin(), contactInterfaces.end());

	std::vector<int> codes(allRegions.size());
	for (size_t r = 0; r < allRegions.size(); r++) {
		BoundaryRegionStore *store = allRegions[r];
		if (store->dimension) {
			for (auto epointer = store->epointers->datatarray().begin(); epointer != store->epointers->datatarray().end(); ++epointer) {
				++store->distribution.code[(int)(*epointer)->code].size;
			}
			for (int code = 0; code < (int)Element::CODE::SIZE; ++code) {
				if (store->distribution.code[code].size) {
					codes[r] |= 1 << code;
				}
			}
		} else {
			store->distribution.code[(int)Element::CODE::POINT1].offset = store->nodeInfo.offset;
		}
	}

	Communication::allReduce(codes.data(), NULL, codes.size(), MPI_INT, MPI_BOR);

	std::vector<esint> sum, offset;
	for (size_t r = 0; r < allRegions.size(); r++) {
		BoundaryRegionStore *store = allRegions[r];
		if (store->dimension) {
			for (size_t i = 0, bitmask = 1; i < store->distribution.code.size(); i++, bitmask = bitmask << 1) {
				if (codes[r] & bitmask) {
					store->distribution.process.size += store->distribution.code[i].size;
					store->distribution.process.next += store->distribution.code[i].size;
					offset.push_back(store->distribution.process.size);
				}
			}
		}
		store->distribution.process.offset = store->distribution.process.size;
		offset.push_back(store->distribution.process.offset);
	}

	sum.resize(offset.size());
	Communication::exscan(sum, offset);

	for (size_t r = 0, j = 0; r < allRegions.size(); r++) {
		BoundaryRegionStore *store = allRegions[r];
		if (store->dimension) {
			for (size_t i = 0, bitmask = 1; i < store->distribution.code.size(); i++, bitmask = bitmask << 1) {
				if (codes[r] & bitmask) {
					store->distribution.code[i].offset = offset[j];
					store->distribution.code[i].totalSize = sum[j++];
				}
			}
		} else {
			store->distribution.code[(int)Element::CODE::POINT1].offset = store->nodeInfo.offset;
			store->distribution.code[(int)Element::CODE::POINT1].totalSize = store->nodeInfo.totalSize;
		}
		store->distribution.process.offset = offset[j];
		store->distribution.process.totalSize = sum[j++];
	}

	profiler::syncend("compute_regions_boudary_distribution");
}

void synchronizeRegionNodes(const NodeStore *nodes, const std::vector<int> &neighbors, std::vector<RegionStore*> &regions)
{
	profiler::syncstart("synchronize_region_nodes");
	std::vector<std::vector<esint> > sBuffer(neighbors.size()), rBuffer(neighbors.size());
	std::vector<size_t> prevsend(neighbors.size());
	for (auto reg = regions.begin(); reg != regions.end(); ++reg) {
		RegionStore *store = *reg;
		for (size_t n = 0; n < prevsend.size(); ++n) {
			prevsend[n] = sBuffer[n].size();
			sBuffer[n].push_back(0);
		}

		esint prev = 0;
		auto ranks = nodes->ranks->begin();
		for (auto n = store->nodes->datatarray().cbegin(); n != store->nodes->datatarray().cend(); prev = *n++) {
			ranks += *n - prev;

			int rindex = 0;
			for (auto r = ranks->begin(); r != ranks->end(); r++) {
				if (*r != info::mpi::rank) {
					while (neighbors[rindex] < *r) { ++rindex; }
					sBuffer[rindex].push_back(nodes->IDs->datatarray()[*n]);
				}
			}
		}

		for (size_t n = 0; n < sBuffer.size(); ++n) {
			std::sort(sBuffer[n].begin() + prevsend[n] + 1, sBuffer[n].end());
		}

		for (size_t n = 0; n < prevsend.size(); ++n) {
			sBuffer[n][prevsend[n]] = sBuffer[n].size() - prevsend[n];
		}
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighbors)) {
		eslog::internalFailure("exchange element region nodes.\n");
	}
	profiler::synccheckpoint("exchange");

	std::fill(prevsend.begin(), prevsend.end(), 0);
	std::vector<size_t> prevrecv(neighbors.size());
	for (auto reg = regions.begin(); reg != regions.end(); ++reg) {
		RegionStore *store = *reg;
		std::vector<esint> nnodes;
		for (size_t n = 0; n < rBuffer.size(); ++n) {
			if (neighbors[n] < info::mpi::rank) {
				for (size_t j = prevrecv[n] + 1, s = prevsend[n] + 1; j < prevrecv[n] + rBuffer[n][prevrecv[n]]; ++j) {
					while (s < prevsend[n] + sBuffer[n][prevsend[n]] && sBuffer[n][s] < rBuffer[n][j]) { ++s; }
					if (s == prevsend[n] + sBuffer[n][prevsend[n]] || sBuffer[n][s] != rBuffer[n][j]) {
						auto it = std::find(nodes->IDs->datatarray().begin(), nodes->IDs->datatarray().begin() + nodes->uniqInfo.nhalo, rBuffer[n][j]);
						nnodes.push_back(it - nodes->IDs->datatarray().begin());
					}
				}
			} else {
				for (size_t j = prevrecv[n] + 1, s = prevsend[n] + 1; j < prevrecv[n] + rBuffer[n][prevrecv[n]]; ++j) {
					while (s < prevsend[n] + sBuffer[n][prevsend[n]] && sBuffer[n][s] < rBuffer[n][j]) { ++s; }
					if (s == prevsend[n] + sBuffer[n][prevsend[n]] || sBuffer[n][s] != rBuffer[n][j]) {
						auto it = std::lower_bound(nodes->IDs->datatarray().begin() + nodes->uniqInfo.nhalo, nodes->IDs->datatarray().end(), rBuffer[n][j]);
						if (it != nodes->IDs->datatarray().end() && *it == rBuffer[n][j]) {
							nnodes.push_back(it - nodes->IDs->datatarray().begin());
						} else {
							auto it = std::find(nodes->IDs->datatarray().begin(), nodes->IDs->datatarray().begin() + nodes->uniqInfo.nhalo, rBuffer[n][j]);
							nnodes.push_back(it - nodes->IDs->datatarray().begin());
						}
					}
				}
			}
			prevrecv[n] += rBuffer[n][prevrecv[n]];
			prevsend[n] += sBuffer[n][prevsend[n]];
		}

		if (nnodes.size()) {
			nnodes.insert(nnodes.end(), store->nodes->datatarray().begin(), store->nodes->datatarray().end());
			utils::sortAndRemoveDuplicates(nnodes);
			delete store->nodes;
			store->nodes = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, nnodes));
		}
	}
	profiler::synccheckpoint("rbuffer");
	profiler::syncend("synchronize_region_nodes");
}

void computeNodeInfo(const NodeStore *nodes, const std::vector<int> &neighbors, std::vector<RegionStore*> &regions)
{
	profiler::syncstart("compute_node_info");
	std::vector<esint> sum, offset;
	for (size_t r = 0; r < regions.size(); ++r) {
		regions[r]->nodeInfo.nhalo = 0;
		for (
				auto n = regions[r]->nodes->datatarray().cbegin();
				n != regions[r]->nodes->datatarray().cend() && *n < nodes->uniqInfo.nhalo;
				++n) {

			++regions[r]->nodeInfo.nhalo;
		}
		regions[r]->nodeInfo.size = regions[r]->nodes->datatarray().size() - regions[r]->nodeInfo.nhalo;
		regions[r]->nodeInfo.offset = regions[r]->nodeInfo.size;
		offset.push_back(regions[r]->nodeInfo.offset);
	}

	sum.resize(offset.size());
	Communication::exscan(sum, offset);
	profiler::synccheckpoint("info");

	for (size_t r = 0, j = 0; r < regions.size(); ++r) {
		regions[r]->nodeInfo.offset = offset[j];
		regions[r]->nodeInfo.totalSize = sum[j++];
	}

	for (size_t r = 0; r < regions.size(); ++r) {
		regions[r]->nodeInfo.position.resize(regions[r]->nodes->datatarray().size());
		std::iota(regions[r]->nodeInfo.position.begin() + regions[r]->nodeInfo.nhalo, regions[r]->nodeInfo.position.end(), regions[r]->nodeInfo.offset);
	}

	std::vector<std::vector<esint> > sBuffer(neighbors.size()), rBuffer(neighbors.size());
	std::vector<size_t> prevsize(neighbors.size()), nsize(neighbors.size());
	std::vector<double> min(3 * regions.size(), std::numeric_limits<double>::max()), max(3 * regions.size(), -std::numeric_limits<double>::max());
	for (size_t r = 0; r < regions.size(); ++r) {
		for (size_t n = 0; n < prevsize.size(); ++n) {
			prevsize[n] = sBuffer[n].size();
			sBuffer[n].push_back(0);
		}

		esint prev = 0, i = 0;
		auto ranks = nodes->ranks->begin();
		for (auto n = regions[r]->nodes->datatarray().cbegin(); n != regions[r]->nodes->datatarray().cend(); prev = *n++, ++i) {
			nodes->coordinates->datatarray()[*n].minmax(min.data(), max.data()); // ALL_ELEMENTS is without nodes
			nodes->coordinates->datatarray()[*n].minmax(min.data() + 3 * r, max.data() + 3 * r);
			ranks += *n - prev;
			if (i < regions[r]->nodeInfo.nhalo) {
				int rindex = 0;
				while (neighbors[rindex] != ranks->front()) { ++rindex; }
				++nsize[rindex];
			} else {
				int rindex = 0;
				for (auto rank = ranks->begin() + 1; rank != ranks->end(); rank++) {
					while (neighbors[rindex] < *rank) { ++rindex; }
					sBuffer[rindex].push_back(regions[r]->nodeInfo.position[i]);
				}
			}
		}

		for (size_t n = 0; n < prevsize.size(); ++n) {
			sBuffer[n][prevsize[n]] = sBuffer[n].size() - prevsize[n];
		}
	}

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		rBuffer[n].resize(nsize[n] + regions.size());
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::receiveLowerKnownSize(sBuffer, rBuffer, neighbors)) {
		eslog::internalFailure("receive global offset of a given element region.\n");
	}
	Communication::allReduce(min, Communication::OP::MIN);
	Communication::allReduce(max, Communication::OP::MAX);
	profiler::synccheckpoint("exchange");

	std::fill(prevsize.begin(), prevsize.end(), 0);
	for (size_t r = 0; r < regions.size(); ++r) {
		regions[r]->nodeInfo.min = Point(min[3 * r], min[3 * r + 1], min[3 * r + 2]);
		regions[r]->nodeInfo.max = Point(max[3 * r], max[3 * r + 1], max[3 * r + 2]);
		for (size_t n = 0, begin = 0; n < rBuffer.size(); ++n) {
			if (neighbors[n] < info::mpi::rank) {
				if (rBuffer[n][prevsize[n]] > 1) {
					std::copy(rBuffer[n].begin() + prevsize[n] + 1, rBuffer[n].begin() + prevsize[n] + rBuffer[n][prevsize[n]], regions[r]->nodeInfo.position.begin() + begin);
					begin += rBuffer[n][prevsize[n]] - 1;
				}
				prevsize[n] += rBuffer[n][prevsize[n]];
			}
		}
	}
	profiler::synccheckpoint("rbuffer");
	profiler::syncend("compute_node_info");
}

void computeRegionsElementNodes(const NodeStore *nodes, const ElementStore *elements, const std::vector<int> &neighbors, std::vector<ElementsRegionStore*> &elementsRegions)
{
	profiler::syncstart("compute_regions_nodes");
	size_t threads = info::env::OMP_NUM_THREADS;
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		std::vector<std::vector<esint> > nodes(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			nodes[t].reserve(elements->nodes->datatarray().size(t));
			auto enodes = elements->nodes->cbegin(t);
			for (auto e = elementsRegions[r]->elements->datatarray().begin(t), prev = e; e != elementsRegions[r]->elements->datatarray().end(t); prev = e++) {
				enodes += *e - *prev;
				nodes[t].insert(nodes[t].end(), enodes->begin(), enodes->end());
			}
			utils::sortAndRemoveDuplicates(nodes[t]);
		}
		utils::mergeThreadedUniqueData(nodes);
		nodes.resize(1);
		nodes.resize(threads);
		serializededata<esint, esint>::balance(1, nodes);

		elementsRegions[r]->nodes = new serializededata<esint, esint>(1, nodes);
	}
	profiler::synccheckpoint("regions_nodes");

	std::vector<RegionStore*> regions(elementsRegions.begin(), elementsRegions.end());
	synchronizeRegionNodes(nodes, neighbors, regions);

	computeNodeInfo(nodes, neighbors, regions);
	profiler::syncend("compute_regions_nodes");
	eslog::checkpointln("MESH: REGIONS NODES COMPUTED");
}

void computeRegionsBoundaryNodes(const std::vector<int> &neighbors, NodeStore *nodes, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces)
{
	profiler::syncstart("compute_regions_boundary_nodes");
	int threads = info::env::OMP_NUM_THREADS;

	std::vector<BoundaryRegionStore*> allRegions;
	allRegions.insert(allRegions.end(), boundaryRegions.begin(), boundaryRegions.end());
	allRegions.insert(allRegions.end(), contactInterfaces.begin(), contactInterfaces.end());

	for (size_t r = 0; r < allRegions.size(); r++) {
		if (allRegions[r]->nodes == NULL) {
			std::vector<std::vector<esint> > nodes(threads);
			nodes[0] = std::vector<esint>(allRegions[r]->elements->datatarray().begin(), allRegions[r]->elements->datatarray().end());
			utils::sortAndRemoveDuplicates(nodes[0]);
			serializededata<esint, esint>::balance(1, nodes);
			allRegions[r]->nodes = new serializededata<esint, esint>(1, nodes);
		}
	}

	profiler::synccheckpoint("compute_nodes");

	std::vector<RegionStore*> regions(allRegions.begin(), allRegions.end());
	synchronizeRegionNodes(nodes, neighbors, regions);

	computeNodeInfo(nodes, neighbors, regions);
	nodes->uniqInfo = allRegions[0]->nodeInfo;
	profiler::syncend("compute_regions_boundary_nodes");
}

void computeBodies(ElementStore *elements, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors)
{
	profiler::syncstart("mesh_bodies_found");

	esint ebegin = elements->distribution.process.offset;
	esint eend = ebegin + elements->distribution.process.size;
	esint nbodies = 0, boffset;
	std::vector<int> body(elements->distribution.process.size, -1);

	{ // DFS for body search
		for (esint e = 0; e < elements->distribution.process.size; ++e) {
			std::vector<esint> stack;
			if (body[e] == -1) {
				stack.push_back(e);
				body[e] = nbodies;
				while (stack.size()) {
					esint current = stack.back();
					stack.pop_back();
					auto dual = elements->faceNeighbors->cbegin() + current;
					for (auto n = dual->begin(); n != dual->end(); ++n) {
						if (*n != -1 && ebegin <= *n && *n < eend) {
							if (body[*n - ebegin] == -1) {
								stack.push_back(*n - ebegin);
								body[*n - ebegin] = nbodies;
							}
						}
					}
				}
				nbodies++;
			}
		}

		// use unique indices across processes
		boffset = nbodies;
		Communication::exscan(boffset);
		for (size_t i = 0; i < body.size(); ++i) {
			body[i] += boffset;
		}
	}

	profiler::synccheckpoint("initial_dfs");

	std::unordered_map<esint, std::unordered_set<esint> > graph; // compact representation (vertex, edges)
	std::unordered_map<esint, esint> holders, labels;
	{ // exchange edges of the compact dual and compute compact representation
		for (esint b = 0; b < nbodies; ++b) {
			graph[b + boffset];
			holders[b + boffset] = info::mpi::rank;
		}
		struct ebody { esint e, b; };
		std::vector<esint> edistribution = Communication::getDistribution(elements->distribution.process.size);
		std::vector<std::vector<ebody> > sBuffer(neighbors.size()), rBuffer(neighbors.size());

		auto e2roffset = [&] (esint e) {
			auto rank = std::lower_bound(edistribution.begin(), edistribution.end(), e + 1) - edistribution.begin() - 1;
			return std::lower_bound(neighbors.begin(), neighbors.end(), rank) - neighbors.begin();
		};

		auto dual = elements->faceNeighbors->cbegin();
		for (esint e = 0; e < elements->distribution.process.size; ++e, ++dual) {
			for (auto n = dual->begin(); n != dual->end(); ++n) {
				if (*n != -1 && (*n < ebegin || eend <= *n)) {
					sBuffer[e2roffset(*n)].push_back(ebody{e + ebegin, body[e]});
				}
			}
		}

		for (size_t n = 0; n < sBuffer.size(); ++n) { // compute compact representation
			std::sort(sBuffer[n].begin(), sBuffer[n].end(), [] (const ebody &i, const ebody &j) { return i.e < j.e; });
			size_t i = 0, last = 0;
			while (i < sBuffer[n].size()) {
				last = i++;
				while (i < sBuffer[n].size() && sBuffer[n][last].b == sBuffer[n][i].b) {
					++i;
				}
				// i points to the first element with different body
				if (sBuffer[n][last].e != sBuffer[n][i - 1].e) {
					sBuffer[n][last] = sBuffer[n][i - 1];
				}
				++last;
			}
			sBuffer[n].resize(last);
		}

		if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighbors)) {
			eslog::internalFailure("cannot exchange bodies indices\n");
		}

		dual = elements->faceNeighbors->cbegin();
		for (esint e = 0; e < elements->distribution.process.size; ++e, ++dual) {
			for (auto n = dual->begin(); n != dual->end(); ++n) {
				if (*n != -1 && (*n < ebegin || eend <= *n)) {
					esint roffset = e2roffset(*n);
					esint nbody = std::lower_bound(rBuffer[roffset].begin(), rBuffer[roffset].end(), *n, [] (const ebody &ebody, esint eindex) { return ebody.e < eindex; })->b;
					graph[body[e]].insert(nbody);
					holders[nbody] = neighbors[roffset];
				}
			}
		}
	}

	profiler::synccheckpoint("graph_synchronization");

	std::unordered_set<int> nn;
	for (auto n = holders.begin(); n != holders.end(); ++n) {
		if (n->second != info::mpi::rank) {
			nn.insert(n->second);
		}
	}

	for (int groups = 1; groups < info::mpi::size; groups = groups << 1) {
		// remove vertices without edges (full body found)
		//////////////////////////////////////////////////
		for (auto v = graph.begin(); v != graph.end(); ) {
			if (v->second.empty()) {
				v = graph.erase(v);
			} else {
				++v;
			}
		}

		// compute the process that receives vertices without edge to other group
		/////////////////////////////////////////////////////////////////////////
		int sink = (info::mpi::rank / groups) % 2 == 0 ? info::mpi::rank + groups : info::mpi::rank - groups;
		if (info::mpi::size <= sink) {
			sink = info::mpi::rank - groups;
		}
		if (0 <= sink && sink < info::mpi::size) {
			nn.insert(sink);
			if (info::mpi::rank + groups < info::mpi::size && info::mpi::size <= info::mpi::rank + 2 * groups) {
				nn.insert(info::mpi::rank + groups);
			}
		} else {
			sink = info::mpi::rank;
		}

		std::vector<int> neighbors(nn.begin(), nn.end());
		std::sort(neighbors.begin(), neighbors.end());
		auto n2i = [&] (int n) {
			return std::lower_bound(neighbors.begin(), neighbors.end(), n) - neighbors.begin();
		};

		// compute target processes
		///////////////////////////
		std::unordered_map<esint, esint> targets = holders;
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				if ((info::mpi::rank / groups) % 2 != (targets[*n] / groups) % 2) {
					if (targets[v->first] == info::mpi::rank) {
						targets[v->first] = targets[*n];
					} else {
						if (std::abs(info::mpi::rank - targets[v->first]) < std::abs(info::mpi::rank - targets[*n])) {
							targets[v->first] = targets[*n];
						}
					}
				}
			}
		}

		// vertices that are without neighbor in the other half are sent to the sink process
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			if (targets[v->first] == info::mpi::rank) {
				targets[v->first] = sink;
			}
		}

		// serialize and exchange target processes
		//////////////////////////////////////////
		std::vector<esint> sTargets;
		std::vector<std::vector<esint> > rTargets(neighbors.size());
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			sTargets.push_back(v->first); sTargets.push_back(targets[v->first]);
		}
		if (!Communication::exchangeUnknownSize(sTargets, rTargets, neighbors)) {
			eslog::internalFailure("cannot exchange compact graph holders.\n");
		}
		for (size_t n = 0; n < rTargets.size(); ++n) {
			for (size_t i = 0; i < rTargets[n].size(); i += 2) {
				targets[rTargets[n][i]] = rTargets[n][i + 1];
			}
		}

		// serialize and exchange graph
		///////////////////////////////
		std::vector<std::vector<esint> > sGraph(neighbors.size()), rGraph(neighbors.size());
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			int t = n2i(targets[v->first]);
			sGraph[t].push_back(v->first);
			sGraph[t].push_back(v->second.size());
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				sGraph[t].push_back(*n);
				if ((targets[*n] / groups) % 2 == (neighbors[t] / groups) % 2) {
					sGraph[t].push_back(targets[*n]);
				} else {
					sGraph[t].push_back(holders[*n]);
				}
			}
		}

		// update holders of exchanged vertices
		///////////////////////////////////////
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				if ((targets[*n] / groups) % 2 == (info::mpi::rank / groups) % 2) {
					holders[*n] = targets[*n];
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sGraph, rGraph, neighbors)) {
			eslog::internalFailure("cannot exchange compact graph.\n");
		}
		for (size_t n = 0; n < rGraph.size(); ++n) {
			for (size_t i = 0; i < rGraph[n].size(); ) {
				esint vertex = rGraph[n][i++];
				holders[vertex] = info::mpi::rank;
				for (esint e = 0; e < rGraph[n][i]; ++e) {
					graph[vertex].insert(rGraph[n][i + 1 + 2 * e]);
					holders[rGraph[n][i + 1 + 2 * e]] = rGraph[n][i + 1 + 2 * e + 1];
				}
				i += 2 * rGraph[n][i] + 1;
			}
		}

		// BFS to compute re-mapping (vertices are renumbered to the lowest number of a neighboring vertex)
		///////////////////////////////////////////////////////////////////////////////////////////////////
		std::unordered_map<esint, esint> remap;
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			esint min = v->first;
			std::unordered_set<esint> marked;
			if (remap.find(v->first) == remap.end()) {
				std::vector<esint> stack({ v->first });
				while (stack.size()) {
					auto e = graph[stack.back()];
					marked.insert(stack.back());
					stack.pop_back();
					for (auto n = e.begin(); n != e.end(); ++n) {
						if (graph.find(*n) != graph.end() && marked.find(*n) == marked.end()) {
							stack.push_back(*n);
							min = std::min(min, *n);
						}
					}
				}
			}
			for (auto m = marked.begin(); m != marked.end(); ++m) {
				remap[*m] = min;
				if (labels.find(*m) != labels.end() || (boffset <= *m && *m < boffset + nbodies)) {
					labels[*m] = min;
					labels[min] = min;
				}
			}
		}

		// exchange re-mapping data in order to synchronize edges
		/////////////////////////////////////////////////////////
		nn.clear();
		for (auto n = holders.begin(); n != holders.end(); ++n) {
			if (n->second != info::mpi::rank) {
				nn.insert(n->second);
			}
		}
		neighbors.assign(nn.begin(), nn.end());
		std::sort(neighbors.begin(), neighbors.end());

		std::vector<esint> sMap;
		std::vector<std::vector<esint> > rMap(neighbors.size());
		for (auto m = remap.begin(); m != remap.end(); ++m) {
			sMap.push_back(m->first); sMap.push_back(m->second);
		}
		if (!Communication::exchangeUnknownSize(sMap, rMap, neighbors)) {
			eslog::internalFailure("cannot exchange re-mapped vertices.\n");
		}
		for (size_t n = 0; n < rMap.size(); ++n) {
			for (size_t i = 0; i < rMap[n].size(); i += 2) {
				remap[rMap[n][i]] = rMap[n][i + 1];
			}
		}

		// update merged graph
		//////////////////////
		std::unordered_map<esint, std::unordered_set<esint> > mgraph;
		std::unordered_map<esint, esint> mholders;
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			mholders[remap[v->first]] = info::mpi::rank;
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				if (graph.find(*n) != graph.end()) {
					if (remap[v->first] != remap[*n]) {
						mgraph[remap[v->first]].insert(remap[*n]);
					}
				} else {
					mgraph[remap[v->first]].insert(remap[*n]);
					mholders[remap[*n]] = holders[*n];
				}
			}
		}

		mgraph.swap(graph);
		mholders.swap(holders);
	}

	for (esint b = boffset; b < boffset + nbodies; ++b) {
		auto it = labels.find(b);
		while (it != labels.end() && it->first != it->second) {
			it = labels.find(it->second);
		}
		labels[b] = it != labels.end() ? it->second : b;
	}

	bodies->size = 0;
	std::vector<esint> ulabels;
	for (esint b = boffset; b < boffset + nbodies; ++b) {
		if (labels[b] == b) {
			ulabels.push_back(b);
			++bodies->size;
		}
	}
	std::sort(ulabels.begin(), ulabels.end());
	bodies->offset = ulabels.size() ? ulabels.front() : 0;
	Communication::allGatherUnknownSize(ulabels);
	bodies->offset = std::lower_bound(ulabels.begin(), ulabels.end(), bodies->offset) - ulabels.begin();

	bodies->totalSize = ulabels.size();
	for (esint b = boffset; b < boffset + nbodies; ++b) {
		labels[b] = std::lower_bound(ulabels.begin(), ulabels.end(), labels[b]) - ulabels.begin();
	}
	for (size_t i = 0; i < body.size(); ++i) {
		body[i] = labels[body[i]];
	}

	if (elements->body == NULL) {
		elements->body = new serializededata<esint, int>(1, tarray<int>(elements->distribution.threads, body));
	} else {
		memcpy(elements->body->datatarray().data(), body.data(), elements->distribution.process.size * sizeof(int));
	}

	int rsize = elements->regions->edataSize();
	std::vector<esint> bodyRegions(bodies->totalSize * rsize);
	for (esint e = 0; e < elements->distribution.process.size; ++e) {
		int b = elements->body->datatarray()[e];
		for (int r = 0; r < rsize; ++r) {
			bodyRegions[rsize * b + r] |= elements->regions->datatarray()[rsize * e + r];
		}
	}

	Communication::allReduce(bodyRegions.data(), NULL, bodyRegions.size(), MPITools::getType<esint>().mpitype, MPI_BOR);

	std::vector<esint> boffsets = { 0 };
	for (esint b = 0; b < bodies->totalSize; ++b) {
		int regions = 0;
		for (int r = 0; r < rsize; ++r) {
			for (size_t bit = 0; bit < 8 * sizeof(esint); ++bit) {
				if (bodyRegions[rsize * b + r] & ((esint)1 << bit)) {
					++regions;
				}
			}
		}
		boffsets.push_back(boffsets.back() + regions); // region ALL_ELEMENTS is not counted
	}

	for (int b = 0; b < bodies->totalSize; ++b) {
		for (int r = 0, rindex = 0; r < rsize; ++r) {
			for (size_t bit = 0; bit < 8 * sizeof(esint) && bit + r * sizeof(esint) < elementsRegions.size(); ++bit, ++rindex) {
				if (bodyRegions[b * rsize + r] & (1 << bit)) {
					elementsRegions[rindex]->bodies.push_back(b);
				}
			}
		}
	}

	std::vector<esint> bcount, esum(bodies->totalSize), fsum(bodies->totalSize);
	for (size_t r = 0; r < elementsRegions.size(); ++r) {
		std::fill(esum.begin(), esum.end(), 0);
		std::fill(fsum.begin(), fsum.end(), 0);
		for (auto e = elementsRegions[r]->elements->datatarray().begin(); e != elementsRegions[r]->elements->datatarray().end(); ++e) {
			int b = elements->body->datatarray()[*e];
			auto neighs = elements->faceNeighbors->begin() + *e;
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n == -1) {
					++fsum[b];
				}
			}
			++esum[b];
		}
		for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
			bcount.push_back(esum[elementsRegions[r]->bodies[b]]);
			bcount.push_back(fsum[elementsRegions[r]->bodies[b]]);
		}
	}

	Communication::allReduce(bcount, Communication::OP::SUM);

	for (size_t r = 0, offset = 0; r < elementsRegions.size(); ++r) {
		for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
			elementsRegions[r]->bodyElements.push_back(bcount[offset++]);
			elementsRegions[r]->bodyFaces.push_back(bcount[offset++]);
		}
	}

	profiler::syncend("mesh_bodies_found");
	eslog::checkpointln("MESH: MESH BODIES FOUND");
}


void triangularizeSurface(SurfaceStore *surface)
{
	if (surface == NULL) {
		return;
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	if (surface->triangles == NULL) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);


		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> ttriangles;
			std::vector<size_t> tintervals;
			if (t == 0) {
				tintervals.push_back(0);
			}

			auto elements = surface->enodes->cbegin(t);
			const auto &epointers = surface->epointers->datatarray().begin();

			for (size_t e = surface->edistribution[t]; e < surface->edistribution[t + 1]; ++e, ++elements) {
				for (auto n = epointers[e]->triangles->datatarray().cbegin(); n != epointers[e]->triangles->datatarray().cend(); ++n) {
					ttriangles.push_back(elements->at(*n));
				}
			}
			tintervals.push_back(ttriangles.size() / 3);

			intervals[t].swap(tintervals);
			triangles[t].swap(ttriangles);
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		surface->tdistribution = intervals[0];
		surface->triangles = new serializededata<esint, esint>(3, triangles);
	}

	eslog::checkpointln("MESH: SURFACE TRIANGULARIZED");
}

void triangularizeBoundary(BoundaryRegionStore *boundary)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	if (boundary->dimension == 2) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);


		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> ttriangles;
			std::vector<size_t> tintervals;
			if (t == 0) {
				tintervals.push_back(0);
			}

			auto elements = boundary->elements->cbegin(t);
			const auto &epointers = boundary->epointers->datatarray().begin();

			for (size_t e = boundary->distribution.threads[t]; e < boundary->distribution.threads[t + 1]; ++e, ++elements) {
				for (auto n = epointers[e]->triangles->datatarray().cbegin(); n != epointers[e]->triangles->datatarray().cend(); ++n) {
					ttriangles.push_back(elements->at(*n));
				}
			}
			tintervals.push_back(ttriangles.size() / 3);

			intervals[t].swap(tintervals);
			triangles[t].swap(ttriangles);
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		boundary->triangles = new serializededata<esint, esint>(3, triangles);
	}

	eslog::checkpointln("MESH: BOUNDARY TRIANGULARIZED");
}

void computeRegionsSurface(ElementStore *elements, NodeStore *nodes, ElementStore *halo, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors)
{
	profiler::syncstart("compute_region_surface");
	size_t threads = info::env::OMP_NUM_THREADS;
	esint eBegin = elements->distribution.process.offset;
	esint eEnd = eBegin + elements->distribution.process.size;

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
		std::vector<std::vector<Element*> > fpointers(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			esint hindex, addFace = 0;
			int rsize = elements->regions->edataSize();
			auto nodes = elements->nodes->cbegin();
			auto neighs = elements->faceNeighbors->cbegin();
			const auto &regions = elements->regions->datatarray();
			const auto &epointers = elements->epointers->datatarray();

			std::vector<esint> fdist, fdata, ecounter((int)Element::CODE::SIZE);
			std::vector<Element*> fpointer;
			if (t == 0) {
				fdist.push_back(0);
			}

			esint prev = 0;
			for (auto e = elementsRegions[r]->elements->datatarray().cbegin(t); e != elementsRegions[r]->elements->datatarray().cend(t); prev = *e++) {
				nodes += *e - prev;
				neighs += *e - prev;
				for (size_t n = 0; n < neighs->size(); ++n) {
					if (neighs->at(n) != -1 && r) {
						if (neighs->at(n) < eBegin || eEnd <= neighs->at(n)) {
							hindex = std::lower_bound(halo->IDs->datatarray().begin(), halo->IDs->datatarray().end(), neighs->at(n)) - halo->IDs->datatarray().begin();
							addFace = memcmp(regions.data() + *e * rsize, halo->regions->datatarray().data() + hindex * rsize, sizeof(esint) * rsize);
						} else {
							addFace = memcmp(regions.data() + *e * rsize, regions.data() + (neighs->at(n) - eBegin) * rsize, sizeof(esint) * rsize);
						}
					} else {
						addFace = neighs->at(n) == -1;
					}
					if (addFace) {
						auto face = epointers[*e]->faces->begin() + n;
						for (auto f = face->begin(); f != face->end(); ++f) {
							fdata.push_back(nodes->at(*f));
						}
						fdist.push_back(fdata.size());
						fpointer.push_back(epointers[*e]->facepointers->datatarray()[n]);
						++ecounter[(int)fpointer.back()->code];
						addFace = 0;
					}
				}
			}

			facesDistribution[t].swap(fdist);
			faces[t].swap(fdata);
			fpointers[t].swap(fpointer);
			ecounters[t].swap(ecounter);
		}

		for (size_t t = 1; t < threads; t++) {
			for (size_t e = 0; e < ecounters[0].size(); e++) {
				ecounters[0][e] += ecounters[t][e];
			}
		}

		serializededata<esint, Element*>::balance(1, fpointers);
		elementsRegions[r]->surface->epointers = new serializededata<esint, Element*>(1, fpointers);
		elementsRegions[r]->surface->ecounters = ecounters[0];

		elementsRegions[r]->surface->edistribution = elementsRegions[r]->surface->epointers->datatarray().distribution();

		if (
				elementsRegions[r]->surface->edistribution.back() &&
				elementsRegions[r]->surface->ecounters[(int)Element::CODE::TRIANGLE3] == (esint)elementsRegions[r]->surface->edistribution.back()) {

			serializededata<esint, esint>::balance(3, faces, &elementsRegions[r]->surface->edistribution);
			elementsRegions[r]->surface->enodes = new serializededata<esint, esint>(3, faces);
			elementsRegions[r]->surface->triangles = elementsRegions[r]->surface->enodes;
			elementsRegions[r]->surface->tdistribution = elementsRegions[r]->surface->edistribution;
		} else {
			utils::threadDistributionToFullDistribution(facesDistribution);
			serializededata<esint, esint>::balance(facesDistribution, faces, &elementsRegions[r]->surface->edistribution);
			elementsRegions[r]->surface->enodes = new serializededata<esint, esint>(facesDistribution, faces);
		}
	}

	profiler::syncend("compute_region_surface");
	eslog::checkpointln("MESH: REGION SURFACE COMPUTED");
}

} // namespace mesh
} // namespace espreso
