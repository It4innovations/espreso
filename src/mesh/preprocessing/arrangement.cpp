
#include "meshpreprocessing.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "mesh/element.h"
#include "mesh/store/store.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/contactinterfacestore.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"
#include "wrappers/metis/w.metis.h"
#include "wrappers/parmetis/w.parmetis.h"

#include <algorithm>
#include <numeric>

namespace espreso {
namespace mesh {


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

	auto localremap = [&] (serializededata<esint, esint>* data) {
		if (data == NULL) {
			return;
		}
		#pragma omp parallel for
		for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
			for (auto e = data->begin(t); e != data->end(t); ++e) {
				for (auto n = e->begin(); n != e->end(); ++n) {
					*n = backpermutation[*n];
				}
			}
		}
	};

	localremap(elements->nodes);

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->elements != NULL) {
			localremap(boundaryRegions[r]->elements);
		}
		if (!StringCompare::caseInsensitiveEq(boundaryRegions[r]->name, "ALL_NODES")) {
			if (boundaryRegions[r]->nodes != NULL) {
				localremap(boundaryRegions[r]->nodes);
				std::sort(boundaryRegions[r]->nodes->datatarray().begin(), boundaryRegions[r]->nodes->datatarray().end());
			}
		}
	}

	profiler::synccheckpoint("remap");
	profiler::syncend("sort_nodes");
	eslog::checkpointln("MESH: NODES ARRANGED");
}

void computeElementIntervals(const DomainStore *domains, ElementStore *elements)
{
	profiler::syncstart("elements_intervals");
	size_t threads = info::env::OMP_NUM_THREADS;
	std::vector<std::vector<esint> > iboundaries(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = domains->distribution[t]; d < domains->distribution[t + 1]; ++d) {
			if (d) {
				iboundaries[t].push_back(domains->elements[d]);
			}
			for (esint e = domains->elements[d] + 1; e < domains->elements[d + 1]; ++e) {
				if (elements->epointers->datatarray()[e]->code != elements->epointers->datatarray()[e - 1]->code) {
					iboundaries[t].push_back(e);
				}
			}
		}
	}
	utils::mergeThreadedUniqueData(iboundaries);
	profiler::synccheckpoint("iboundaries");

	elements->eintervals.clear();
	elements->eintervals.push_back(ElementsInterval(0, 0));
	elements->eintervals.back().domain = domains->offset;
	elements->eintervals.back().code = static_cast<int>(elements->epointers->datatarray().front()->code);
	elements->eintervalsDistribution.push_back(0);
	for (size_t i = 0; i < iboundaries[0].size(); i++) {
		elements->eintervals.back().end = iboundaries[0][i];
		elements->eintervals.push_back(ElementsInterval(iboundaries[0][i], iboundaries[0][i]));
		const std::vector<esint> &edist = domains->elements;
		elements->eintervals.back().domain = std::lower_bound(edist.begin(), edist.end(), elements->eintervals.back().begin + 1) - edist.begin() - 1 + domains->offset;
		elements->eintervals.back().code = static_cast<int>(elements->epointers->datatarray()[elements->eintervals.back().begin]->code);
		if ((elements->eintervals.end() - 1)->domain != (elements->eintervals.end() - 2)->domain) {
			elements->eintervalsDistribution.push_back(elements->eintervals.size() - 1);
		}
	}
	elements->eintervals.back().end = elements->distribution.process.size;
	elements->eintervalsDistribution.push_back(elements->eintervals.size());
	profiler::syncend("elements_intervals");
	eslog::checkpointln("MESH: ELEMENTS INTERVALS COMPUTED");
}

void computeElementDistribution(const std::vector<ElementsInterval> &eintervals, ElementsDistributionInfo &distribution)
{
	profiler::syncstart("compute_element_distribution");
	distribution.clear();
	for (size_t i = 0; i < eintervals.size(); ++i) {
		distribution.process.size += eintervals[i].end - eintervals[i].begin;
		distribution.code[eintervals[i].code].size += eintervals[i].end - eintervals[i].begin;
	}

	std::vector<esint> sum, offset;
	for (size_t i = 0; i < distribution.code.size(); ++i) {
		offset.push_back(distribution.code[i].size);
	}
	sum.resize(offset.size());
	Communication::exscan(sum, offset);
	for (size_t i = 0; i < distribution.code.size(); ++i) {
		distribution.code[i].offset = offset[i];
		distribution.code[i].totalSize = sum[i];
	}

	distribution.process.offset = distribution.process.size;
	distribution.process.totalSize = Communication::exscan(distribution.process.offset);
	distribution.process.next = distribution.process.offset + distribution.process.size;
	profiler::syncend("compute_element_distribution");
}

void computeRegionsElementIntervals(const ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions)
{
	profiler::syncstart("compute_regions_elements_intervals");

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		const auto &relements = elementsRegions[r]->elements->datatarray();
		elementsRegions[r]->eintervals = elements->eintervals;
		for (size_t i = 0; i < elementsRegions[r]->eintervals.size(); ++i) {
			elementsRegions[r]->eintervals[i].begin = std::lower_bound(relements.begin(), relements.end(), elementsRegions[r]->eintervals[i].begin) - relements.begin();
			elementsRegions[r]->eintervals[i].end = std::lower_bound(relements.begin(), relements.end(), elementsRegions[r]->eintervals[i].end) - relements.begin();
		}
	}
	profiler::syncend("compute_regions_elements_intervals");
	eslog::checkpointln("MESH: ELEMENTS REGIONS INTERVALS COMPUTED");
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

void computeRegionsElementDistribution(const ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions)
{
	profiler::syncstart("regions element distribution");
	std::vector<esint> sum, offset;
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		elementsRegions[r]->distribution.clear();
		for (size_t i = 0; i < elementsRegions[r]->eintervals.size(); ++i) {
			elementsRegions[r]->distribution.code[elementsRegions[r]->eintervals[i].code].size += elementsRegions[r]->eintervals[i].end - elementsRegions[r]->eintervals[i].begin;
		}
		for (size_t i = 0; i < elements->distribution.code.size(); i++) {
			if (elements->distribution.code[i].totalSize) {
				elementsRegions[r]->distribution.process.size += elementsRegions[r]->distribution.code[i].size;
				elementsRegions[r]->distribution.process.next += elementsRegions[r]->distribution.code[i].size;
				offset.push_back(elementsRegions[r]->distribution.code[i].size);
			}
		}
		elementsRegions[r]->distribution.process.offset = elementsRegions[r]->distribution.process.size;
		offset.push_back(elementsRegions[r]->distribution.process.offset);
	}

	sum.resize(offset.size());
	Communication::exscan(sum, offset);

	for (size_t r = 0, j = 0; r < elementsRegions.size(); r++) {
		for (size_t i = 0; i < elements->distribution.code.size(); i++) {
			if (elements->distribution.code[i].totalSize) {
				elementsRegions[r]->distribution.code[i].offset = offset[j];
				elementsRegions[r]->distribution.code[i].totalSize = sum[j++];
			}
		}
		elementsRegions[r]->distribution.process.offset = offset[j];
		elementsRegions[r]->distribution.process.next += offset[j];
		elementsRegions[r]->distribution.process.totalSize = sum[j++];
	}
	profiler::syncend("regions element distribution");
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

void computeRegionsBoundaryParents(const NodeStore *nodes, const ElementStore *elements, const DomainStore *domains, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces)
{
	profiler::syncstart("arrange_boudary_regions");
	int threads = info::env::OMP_NUM_THREADS;

	std::vector<BoundaryRegionStore*> allRegions;
	allRegions.insert(allRegions.end(), boundaryRegions.begin(), boundaryRegions.end());
	allRegions.insert(allRegions.end(), contactInterfaces.begin(), contactInterfaces.end());

	esint eoffset = elements->distribution.process.offset;
	for (size_t r = 0; r < allRegions.size(); r++) {
		BoundaryRegionStore *store = allRegions[r];
		if (store->dimension == 0) {
			store->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, 0));
			store->epointers = new serializededata<esint, Element*>(1, tarray<Element*>(threads, 0));
		} else {
			std::vector<size_t> distribution = tarray<size_t>::distribute(threads, store->elements->structures());
			const std::vector<esint> &eDomainDistribution = domains->elements;
			std::vector<esint> emembership(distribution.back()), edomain(distribution.back());

			#pragma omp parallel for
			for (int t = 0; t < threads; t++) {
				auto enodes = store->elements->cbegin() + distribution[t];
				std::vector<esint> nlinks;
				size_t counter;
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					nlinks.clear();
					for (auto n = enodes->begin(); n != enodes->end(); ++n) {
						auto links = nodes->elements->cbegin() + *n;
						nlinks.insert(nlinks.end(), links->begin(), links->end());
					}
					std::sort(nlinks.begin(), nlinks.end());
					counter = 1;
					for (size_t i = 1; i < nlinks.size(); ++i) {
						if (nlinks[i - 1] == nlinks[i]) {
							++counter;
							if (counter == enodes->size() && eoffset <= nlinks[i]) {
								emembership[e] = nlinks[i];
								edomain[e] = std::lower_bound(eDomainDistribution.begin(), eDomainDistribution.end(), nlinks[i] + 1 - eoffset) - eDomainDistribution.begin() - 1;
								break;
							}
						} else {
							counter = 1;
						}
					}
				}
			}

			std::vector<esint> permutation(store->elements->structures());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
				if (edomain[i] == edomain[j]) {
					if (store->epointers->datatarray()[i] == store->epointers->datatarray()[j]) {
						return emembership[i] < emembership[j];
					}
					return store->epointers->datatarray()[i] < store->epointers->datatarray()[j];
				}
				return edomain[i] < edomain[j];
			});

			std::vector<std::vector<esint> > ememberdata(threads);

			#pragma omp parallel for
			for (int t = 0; t < threads; t++) {
				std::vector<esint> tdata;

				auto nodes = store->elements->begin() + distribution[t];
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++nodes) {
					esint eindex = emembership[e] - eoffset;
					auto epointer = elements->epointers->datatarray()[eindex];
					auto enodes = elements->nodes->begin() + eindex;
					tdata.push_back(eindex);
					if (store->dimension == 1) {
						tdata.push_back(epointer->getIndex(*enodes, epointer->edges, epointer->edgepointers, *nodes));
					}
					if (store->dimension == 2) {
						tdata.push_back(epointer->getIndex(*enodes, epointer->faces, epointer->facepointers, *nodes));
					}
					if (tdata.back() < 0) {
						eslog::internalFailure("cannot find sub-element index.\n");
					}
				}

				ememberdata[t].swap(tdata);
			}

			store->emembership = new serializededata<esint, esint>(2, ememberdata);

			std::vector<size_t> edistribution;
			for (auto i = domains->elements.begin(); i != domains->elements.end(); ++i) {
				auto it = std::lower_bound(permutation.begin(), permutation.end(), *i, [&] (esint i, esint d) { return emembership[i] - eoffset < d; });
				edistribution.push_back(it - permutation.begin());
			}

			std::vector<size_t> tdistribution;
			for (size_t t = 0; t < domains->distribution.size(); t++) {
				tdistribution.push_back(edistribution[domains->distribution[t]]);
			}

			store->permute(permutation, tdistribution);

			std::vector<std::vector<esint> > iboundaries(threads);

			#pragma omp parallel for
			for (int t = 0; t < threads; t++) {
				for (size_t d = domains->distribution[t]; d < domains->distribution[t + 1]; d++) {
					iboundaries[t].push_back(edistribution[d]);
					for (size_t e = edistribution[d] + 1; e < edistribution[d + 1]; ++e) {
						if (store->epointers->datatarray()[e]->code != store->epointers->datatarray()[e - 1]->code) {
							iboundaries[t].push_back(e);
						}
					}
				}
			}
			iboundaries.back().push_back(edistribution.back());
			utils::mergeThreadedUniqueData(iboundaries);

			store->eintervalsDistribution.push_back(0);
			esint lastDomain = 0;
			for (size_t i = 0; i < iboundaries[0].size() - 1; i++) {
				store->eintervals.push_back(ElementsInterval(iboundaries[0][i], iboundaries[0][i + 1]));
				store->eintervals.back().code = static_cast<int>(store->epointers->datatarray()[iboundaries[0][i]]->code);
				store->eintervals.back().domain = edomain[permutation[iboundaries[0][i]]];
				if (store->eintervals.back().domain != lastDomain) {
					store->eintervalsDistribution.insert(
							store->eintervalsDistribution.end(),
							store->eintervals.back().domain - lastDomain,
							store->eintervals.size() - 1);
				}
				lastDomain = store->eintervals.back().domain;
			}
			store->eintervalsDistribution.insert(
					store->eintervalsDistribution.end(),
					domains->size - lastDomain,
					store->eintervals.size());
		}
	}

	profiler::syncend("arrange_boudary_regions");
	eslog::checkpointln("MESH: BOUNDARY REGIONS ARRANGED");
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
			for (size_t i = 0; i < store->eintervals.size(); ++i) {
				store->distribution.code[store->eintervals[i].code].size += store->eintervals[i].end - store->eintervals[i].begin;
				codes[r] |= 1 << store->eintervals[i].code;
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

void fillRegionMask(const ElementsDistributionInfo &distribution, const std::vector<ElementsRegionStore*> &elementsRegions, serializededata<esint, esint>* &mask)
{
	if (mask != NULL) {
		return;
	}
	profiler::syncstart("fill_region_mask");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > eregions(threads);

	// regions are transfered via mask
	int regionsBitMaskSize = bitMastSize(elementsRegions.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		eregions[t].resize(regionsBitMaskSize * (distribution.threads[t + 1] - distribution.threads[t]));
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));

			const auto &elements = elementsRegions[r]->elements->datatarray();
			auto begin = std::lower_bound(elements.begin(), elements.end(), distribution.threads[t]);
			auto end = std::lower_bound(elements.begin(), elements.end(), distribution.threads[t + 1]);
			for (auto i = begin; i != end; ++i) {
				eregions[t][(*i - distribution.threads[t]) * regionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}
	mask = new serializededata<esint, esint>(regionsBitMaskSize, eregions);

	profiler::syncend("fill_region_mask");
	eslog::checkpointln("MESH: REGION MASK FILLED");
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

}
}

