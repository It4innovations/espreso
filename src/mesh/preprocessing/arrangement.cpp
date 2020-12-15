
#include "meshpreprocessing.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "mesh/element.h"
#include "mesh/store/store.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"
#include "wrappers/metis/w.metis.h"
#include "wrappers/parmetis/w.parmetis.h"

#include <algorithm>
#include <numeric>

namespace espreso {
namespace mesh {

void arrangeNodes()
{
	profiler::syncstart("arange_nodes");
	profiler::syncparam("size", info::mesh->nodes->size);
	std::vector<esint> permutation(info::mesh->nodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		int irank = *(info::mesh->nodes->ranks->begin() + i)->begin();
		int jrank = *(info::mesh->nodes->ranks->begin() + j)->begin();
		if (irank == jrank) {
			return info::mesh->nodes->IDs->datatarray()[i] < info::mesh->nodes->IDs->datatarray()[j];
		}
		return irank < jrank;
	});

	info::mesh->nodes->permute(permutation);
	profiler::synccheckpoint("permute");

	// nhalo is used by other routines (before arrange boundary elements)
	info::mesh->nodes->uniqInfo.nhalo = 0;
	auto ranks = info::mesh->nodes->ranks->begin();
	while (ranks != info::mesh->nodes->ranks->end() && *ranks->begin() < info::mpi::rank) {
		++ranks;
		++info::mesh->nodes->uniqInfo.nhalo;
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

	localremap(info::mesh->elements->procNodes);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->procNodes != NULL) {
			localremap(info::mesh->boundaryRegions[r]->procNodes);
		}
		if (!StringCompare::caseInsensitiveEq(info::mesh->boundaryRegions[r]->name, "ALL_NODES")) {
			if (info::mesh->boundaryRegions[r]->nodes != NULL) {
				localremap(info::mesh->boundaryRegions[r]->nodes);
				std::sort(info::mesh->boundaryRegions[r]->nodes->datatarray().begin(), info::mesh->boundaryRegions[r]->nodes->datatarray().end());
			}
		}
	}

	profiler::synccheckpoint("remap");
	profiler::syncend("arange_nodes");
	eslog::checkpointln("MESH: NODES ARRANGED");
}

void arrangeElements()
{
	profiler::syncstart("arrange_elements");
	std::vector<esint> permutation(info::mesh->elements->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	arrangeElementsPermutation(permutation);
	permuteElements(permutation, info::mesh->elements->distribution);

	profiler::syncend("arrange_elements");
	eslog::checkpointln("MESH: ELEMENTS ARRANGED");
}

void arrangeElementsPermutation(std::vector<esint> &permutation)
{
	profiler::syncstart("permute_elements");
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; ++d) {
			std::sort(
					permutation.begin() + info::mesh->elements->elementsDistribution[d],
					permutation.begin() + info::mesh->elements->elementsDistribution[d + 1],
					[&] (esint i, esint j) {
				if (info::mesh->elements->epointers->datatarray()[i]->code != info::mesh->elements->epointers->datatarray()[j]->code) {
					return info::mesh->elements->epointers->datatarray()[i]->code < info::mesh->elements->epointers->datatarray()[j]->code;
				}
				return i < j;
			});
		}
	}

	profiler::synccheckpoint("sort_etypes");

	std::vector<std::vector<esint> > iboundaries(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; ++d) {
			if (d) {
				iboundaries[t].push_back(info::mesh->elements->elementsDistribution[d]);
			}
			for (esint e = info::mesh->elements->elementsDistribution[d] + 1; e < info::mesh->elements->elementsDistribution[d + 1]; ++e) {
				if (info::mesh->elements->epointers->datatarray()[permutation[e]]->code != info::mesh->elements->epointers->datatarray()[permutation[e - 1]]->code) {
					iboundaries[t].push_back(e);
				}
			}
		}
	}
	utils::mergeThreadedUniqueData(iboundaries);
	profiler::synccheckpoint("iboundaries");

	info::mesh->elements->eintervals.push_back(ElementsInterval(0, 0));
	info::mesh->elements->eintervals.back().domain = info::mesh->elements->firstDomain;
	info::mesh->elements->eintervals.back().code = static_cast<int>(info::mesh->elements->epointers->datatarray()[permutation[0]]->code);
	info::mesh->elements->eintervalsDistribution.push_back(0);
	for (size_t i = 0; i < iboundaries[0].size(); i++) {
		info::mesh->elements->eintervals.back().end = iboundaries[0][i];
		info::mesh->elements->eintervals.push_back(ElementsInterval(iboundaries[0][i], iboundaries[0][i]));
		const std::vector<esint> &edist = info::mesh->elements->elementsDistribution;
		info::mesh->elements->eintervals.back().domain = std::lower_bound(edist.begin(), edist.end(), info::mesh->elements->eintervals.back().begin + 1) - edist.begin() - 1 + info::mesh->elements->firstDomain;
		info::mesh->elements->eintervals.back().code = static_cast<int>(info::mesh->elements->epointers->datatarray()[permutation[info::mesh->elements->eintervals.back().begin]]->code);
		if ((info::mesh->elements->eintervals.end() - 1)->domain != (info::mesh->elements->eintervals.end() - 2)->domain) {
			info::mesh->elements->eintervalsDistribution.push_back(info::mesh->elements->eintervals.size() - 1);
		}
	}
	info::mesh->elements->eintervals.back().end = info::mesh->elements->size;
	info::mesh->elements->eintervalsDistribution.push_back(info::mesh->elements->eintervals.size());
	profiler::synccheckpoint("eintervals");

	int elementstypes = static_cast<int>(Element::CODE::SIZE);
	if (elementstypes > 32) {
		eslog::error("ESPRESO internal error: increase elements-types synchronization buffer.\n");
	}

	int codes = 0;
	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		codes |= 1 << info::mesh->elements->eintervals[i].code;
		info::mesh->elements->ecounters[info::mesh->elements->eintervals[i].code] += info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
	}

	int allcodes = 0;
	Communication::allReduce(&codes, &allcodes, 1, MPI_INT, MPI_BOR);

	std::vector<esint> sum, offset;
	for (int i = 0, bitmask = 1; i < elementstypes; i++, bitmask = bitmask << 1) {
		if (allcodes & bitmask) {
			offset.push_back(info::mesh->elements->ecounters[i]);
		}
	}
	sum.resize(offset.size());
	Communication::exscan(sum, offset);
	for (int i = 0, bitmask = 1, j = 0; i < elementstypes; i++, bitmask = bitmask << 1) {
		if (allcodes & bitmask) {
			info::mesh->elements->ecounters[i] = sum[j++];
		}
	}

	std::vector<esint> procdist = Communication::getDistribution(info::mesh->elements->size);
	info::mesh->elements->offset = procdist[info::mpi::rank];
	info::mesh->elements->totalSize = procdist.back();

	profiler::synccheckpoint("elements_statistics");
	profiler::syncend("permute_elements");
	eslog::checkpointln("MESH: ELEMENTS PERMUTATION ARRANGED");
}

void arrangeElementsRegions()
{
	profiler::syncstart("arrange_elements_regions");
	if (info::mesh->elements->regions == NULL) {
		fillRegionMask();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
		const auto &elements = info::mesh->elementsRegions[r]->elements->datatarray();

		info::mesh->elementsRegions[r]->eintervals = info::mesh->elements->eintervals;
		for (size_t i = 0; i < info::mesh->elementsRegions[r]->eintervals.size(); ++i) {
			info::mesh->elementsRegions[r]->eintervals[i].begin = std::lower_bound(elements.begin(), elements.end(), info::mesh->elementsRegions[r]->eintervals[i].begin) - elements.begin();
			info::mesh->elementsRegions[r]->eintervals[i].end = std::lower_bound(elements.begin(), elements.end(), info::mesh->elementsRegions[r]->eintervals[i].end) - elements.begin();
		}
	}

	profiler::synccheckpoint("eintervals");

	for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
		std::vector<std::vector<esint> > unique(threads);
		info::mesh->elementsRegions[r]->ueintervals = info::mesh->elementsRegions[r]->eintervals;

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			const auto &regions = info::mesh->elements->regions->datatarray();
			int maskSize = info::mesh->elements->regionMaskSize;
			esint maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));
			std::vector<esint> mask(maskSize);
			mask[maskOffset] = bit;

			for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
				for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; i++) {
					size_t usize = unique[t].size();
					for (esint e = info::mesh->elementsRegions[r]->eintervals[i].begin; e < info::mesh->elementsRegions[r]->eintervals[i].end; ++e) {
						if (memcmp(regions.data() + e * maskSize, mask.data(), sizeof(esint) * maskSize) == 0) {
							unique[t].push_back(e);
						}
					}
					info::mesh->elementsRegions[r]->ueintervals[i].begin = 0;
					info::mesh->elementsRegions[r]->ueintervals[i].end = unique[t].size() - usize;
				}
			}
		}

		for (size_t i = 1; i < info::mesh->elementsRegions[r]->ueintervals.size(); ++i) {
			info::mesh->elementsRegions[r]->ueintervals[i].begin += info::mesh->elementsRegions[r]->ueintervals[i - 1].end;
			info::mesh->elementsRegions[r]->ueintervals[i].end   += info::mesh->elementsRegions[r]->ueintervals[i - 1].end;
		}

		if (info::mesh->elementsRegions[r]->eintervals == info::mesh->elementsRegions[r]->ueintervals) {
			info::mesh->elementsRegions[r]->uniqueElements = info::mesh->elementsRegions[r]->elements;
		} else {
			info::mesh->elementsRegions[r]->uniqueElements = new serializededata<esint, esint>(1, unique);
		}

		auto computeCountersAndNodes = [&] (const std::vector<ElementsInterval> &eintervals, const tarray<esint> &elements) {
			for (size_t i = 0; i < eintervals.size(); ++i) {
				info::mesh->elementsRegions[r]->ecounters[eintervals[i].code] += eintervals[i].end - eintervals[i].begin;
			}

			std::vector<std::vector<esint> > nodes(threads);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
					for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; i++) {
						if (eintervals[i].end - eintervals[i].begin > 0) {
							if (eintervals[i].end - eintervals[i].begin == info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) {
								nodes[t].insert(
										nodes[t].end(),
										(info::mesh->elements->procNodes->cbegin() + info::mesh->elements->eintervals[i].begin)->begin(),
										(info::mesh->elements->procNodes->cbegin() + info::mesh->elements->eintervals[i].end)->begin());
							} else {
								auto enodes = info::mesh->elements->procNodes->cbegin() + info::mesh->elements->eintervals[i].begin;
								esint prev = info::mesh->elements->eintervals[i].begin;
								for (esint e = eintervals[i].begin; e < eintervals[i].end; prev = elements[e++]) {
									enodes += elements[e] - prev;
									nodes[t].insert(nodes[t].end(), enodes->begin(), enodes->end());
								}
							}
						}
					}
				}
				utils::sortAndRemoveDuplicates(nodes[t]);
			}
			utils::mergeThreadedUniqueData(nodes);
			nodes.resize(1);
			nodes.resize(threads);
			serializededata<esint, esint>::balance(1, nodes);

			info::mesh->elementsRegions[r]->nodes = new serializededata<esint, esint>(1, nodes);
		};

		if (StringCompare::caseInsensitiveEq(info::mesh->elementsRegions[r]->name, "ALL_ELEMENTS")) {
			computeCountersAndNodes(info::mesh->elementsRegions[r]->ueintervals, info::mesh->elementsRegions[r]->uniqueElements->datatarray());
		} else {
			computeCountersAndNodes(info::mesh->elementsRegions[r]->eintervals, info::mesh->elementsRegions[r]->elements->datatarray());
		}
	}
	profiler::synccheckpoint("regions_nodes");

	std::vector<RegionStore*> regions(info::mesh->elementsRegions.begin(), info::mesh->elementsRegions.end());
	synchronizeRegionNodes(regions);

	computeNodeInfo(regions);
	profiler::synccheckpoint("node_info");

	std::vector<esint> sum, offset;
	for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
		for (size_t i = 0; i < info::mesh->elements->ecounters.size(); i++) {
			if (info::mesh->elements->ecounters[i]) {
				info::mesh->elementsRegions[r]->size += info::mesh->elementsRegions[r]->ecounters[i];
				offset.push_back(info::mesh->elementsRegions[r]->ecounters[i]);
			}
		}
		info::mesh->elementsRegions[r]->offset = info::mesh->elementsRegions[r]->size;
		offset.push_back(info::mesh->elementsRegions[r]->offset);
	}

	sum.resize(offset.size());
	Communication::exscan(sum, offset);

	for (size_t r = 0, j = 0; r < info::mesh->elementsRegions.size(); r++) {
		for (size_t i = 0; i < info::mesh->elements->ecounters.size(); i++) {
			if (info::mesh->elements->ecounters[i]) {
				info::mesh->elementsRegions[r]->eoffsets[i] = offset[j];
				info::mesh->elementsRegions[r]->ecounters[i] = sum[j++];
			}
		}
		info::mesh->elementsRegions[r]->offset = offset[j];
		info::mesh->elementsRegions[r]->totalsize = sum[j++];
	}

	if (info::mesh->eregion("ALL_ELEMENTS")->nodeInfo.totalSize) {
		ElementsRegionStore* nameless = info::mesh->eregion("ALL_ELEMENTS");
		info::mesh->elementsRegions.push_back(new ElementsRegionStore("NAMELESS_ELEMENT_SET"));
		info::mesh->elementsRegions.back()->size = nameless->size;
		info::mesh->elementsRegions.back()->offset = nameless->offset;
		info::mesh->elementsRegions.back()->totalsize = nameless->totalsize;
		info::mesh->elementsRegions.back()->eoffsets = nameless->eoffsets;
		info::mesh->elementsRegions.back()->ecounters = nameless->ecounters;
		info::mesh->elementsRegions.back()->eintervals = nameless->ueintervals;
		info::mesh->elementsRegions.back()->elements = new serializededata<esint, esint>(*nameless->uniqueElements);
		info::mesh->elementsRegions.back()->uniqueElements = info::mesh->elementsRegions.back()->elements;
		info::mesh->elementsRegions.back()->nodes = new serializededata<esint, esint>(*nameless->nodes);
		info::mesh->elementsRegions.back()->nodeInfo.nhalo = nameless->nodeInfo.nhalo;
		info::mesh->elementsRegions.back()->nodeInfo.offset = nameless->nodeInfo.offset;
		info::mesh->elementsRegions.back()->nodeInfo.size = nameless->nodeInfo.size;
		info::mesh->elementsRegions.back()->nodeInfo.totalSize = nameless->nodeInfo.totalSize;
		info::mesh->elementsRegions.back()->nodeInfo.position = nameless->nodeInfo.position;
	}

	profiler::synccheckpoint("element_info");
	profiler::syncend("arrange_elements_regions");
	eslog::checkpointln("MESH: ELEMENTS REGIONS ARRANGED");
}

void arrangeBoundaryRegions()
{
	profiler::syncstart("arrange_boudary_regions");
	int threads = info::env::OMP_NUM_THREADS;

	esint eoffset = info::mesh->elements->offset;
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->nodes == NULL) {
			std::vector<std::vector<esint> > nodes(threads);
			nodes[0] = std::vector<esint>(info::mesh->boundaryRegions[r]->procNodes->datatarray().begin(), info::mesh->boundaryRegions[r]->procNodes->datatarray().end());
			utils::sortAndRemoveDuplicates(nodes[0]);
			serializededata<esint, esint>::balance(1, nodes);
			info::mesh->boundaryRegions[r]->nodes = new serializededata<esint, esint>(1, nodes);
		}
	}

	profiler::synccheckpoint("compute_nodes");

	std::vector<RegionStore*> regions(info::mesh->boundaryRegions.begin(), info::mesh->boundaryRegions.end());
	synchronizeRegionNodes(regions);

	computeNodeInfo(regions);
	info::mesh->nodes->uniqInfo = info::mesh->boundaryRegions[0]->nodeInfo;
	profiler::synccheckpoint("node_info");

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *store = info::mesh->boundaryRegions[r];
		if (store->dimension == 0) {
			store->procNodes = new serializededata<esint, esint>(1, tarray<esint>(threads, 0));
			store->epointers = new serializededata<esint, Element*>(1, tarray<Element*>(threads, 0));
		} else {
			std::vector<size_t> distribution = tarray<size_t>::distribute(threads, store->procNodes->structures());
			std::vector<esint> &eDomainDistribution = info::mesh->elements->elementsDistribution;
			std::vector<esint> emembership(distribution.back()), edomain(distribution.back());

			#pragma omp parallel for
			for (int t = 0; t < threads; t++) {
				auto enodes = store->procNodes->cbegin() + distribution[t];
				std::vector<esint> nlinks;
				size_t counter;
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					nlinks.clear();
					for (auto n = enodes->begin(); n != enodes->end(); ++n) {
						auto links = info::mesh->nodes->elements->cbegin() + *n;
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

			std::vector<esint> permutation(store->procNodes->structures());
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

				auto nodes = store->procNodes->begin() + distribution[t];
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++nodes) {
					esint eindex = emembership[e] - eoffset;
					auto epointer = info::mesh->elements->epointers->datatarray()[eindex];
					auto enodes = info::mesh->elements->procNodes->begin() + eindex;
					tdata.push_back(eindex);
					if (store->dimension == 1) {
						tdata.push_back(epointer->getIndex(*enodes, epointer->edges, epointer->edgepointers, *nodes));
					}
					if (store->dimension == 2) {
						tdata.push_back(epointer->getIndex(*enodes, epointer->faces, epointer->facepointers, *nodes));
					}
					if (tdata.back() < 0) {
						eslog::error("ESPRESO internal error: cannot find sub-element index.\n");
					}
				}

				ememberdata[t].swap(tdata);
			}

			store->emembership = new serializededata<esint, esint>(2, ememberdata);

			std::vector<size_t> edistribution;
			for (auto i = info::mesh->elements->elementsDistribution.begin(); i != info::mesh->elements->elementsDistribution.end(); ++i) {
				auto it = std::lower_bound(permutation.begin(), permutation.end(), *i, [&] (esint i, esint d) { return emembership[i] - eoffset < d; });
				edistribution.push_back(it - permutation.begin());
			}

			std::vector<size_t> tdistribution;
			for (size_t t = 0; t < info::mesh->elements->domainDistribution.size(); t++) {
				tdistribution.push_back(edistribution[info::mesh->elements->domainDistribution[t]]);
			}

			store->permute(permutation, tdistribution);

			std::vector<std::vector<esint> > iboundaries(threads);

			#pragma omp parallel for
			for (int t = 0; t < threads; t++) {
				for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
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
					info::mesh->elements->ndomains - lastDomain,
					store->eintervals.size());
		}
	}

	profiler::synccheckpoint("emembership");

	std::vector<int> codes(info::mesh->boundaryRegions.size()), gcodes(info::mesh->boundaryRegions.size());
	esint bsize = 0;
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *store = info::mesh->boundaryRegions[r];
		if (info::mesh->boundaryRegions[r]->dimension) {
			bsize += info::mesh->boundaryRegions[r]->epointers->datatarray().size();

			for (size_t i = 0; i < store->eintervals.size(); ++i) {
				store->ecounters[store->eintervals[i].code] += store->eintervals[i].end - store->eintervals[i].begin;
				codes[r] |= 1 << store->eintervals[i].code;
			}
		} else {
			store->ecounters[(int)Element::CODE::POINT1] = store->nodeInfo.offset;
		}
	}

	Communication::allReduce(codes.data(), gcodes.data(), codes.size(), MPI_INT, MPI_BOR);

	std::vector<esint> sum, offset;
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *store = info::mesh->boundaryRegions[r];
		if (info::mesh->boundaryRegions[r]->dimension) {
			for (size_t i = 0, bitmask = 1; i < info::mesh->elements->ecounters.size(); i++, bitmask = bitmask << 1) {
				if (gcodes[r] & bitmask) {
					store->size += store->ecounters[i];
					offset.push_back(store->ecounters[i]);
				}
			}
		}
		store->offset = store->size;
		offset.push_back(store->offset);
	}

	sum.resize(offset.size());
	Communication::exscan(sum, offset);

	for (size_t r = 0, j = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *store = info::mesh->boundaryRegions[r];
		if (info::mesh->boundaryRegions[r]->dimension) {
			for (size_t i = 0, bitmask = 1; i < info::mesh->elements->ecounters.size(); i++, bitmask = bitmask << 1) {
				if (gcodes[r] & bitmask) {
					store->eoffsets[i] = offset[j];
					store->ecounters[i] = sum[j++];
				}
			}
		} else {
			store->eoffsets[(int)Element::CODE::POINT1] = store->nodeInfo.offset;
			store->ecounters[(int)Element::CODE::POINT1] = store->nodeInfo.totalSize;
		}
		store->offset = offset[j];
		store->totalsize = sum[j++];
	}

	profiler::synccheckpoint("boundary_stats");
	profiler::syncend("arrange_boudary_regions");
	eslog::checkpointln("MESH: BOUNDARY REGIONS ARRANGED");
}

void fillRegionMask()
{
	profiler::syncstart("fill_region_mask");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > eregions(threads);

	// regions are transfered via mask
	int regionsBitMaskSize = info::mesh->elementsRegions.size() / (8 * sizeof(esint)) + (info::mesh->elementsRegions.size() % (8 * sizeof(esint)) ? 1 : 0);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		eregions[t].resize(regionsBitMaskSize * (info::mesh->elements->distribution[t + 1] - info::mesh->elements->distribution[t]));
		for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));

			const auto &elements = info::mesh->elementsRegions[r]->elements->datatarray();
			auto begin = std::lower_bound(elements.begin(), elements.end(), info::mesh->elements->distribution[t]);
			auto end = std::lower_bound(elements.begin(), elements.end(), info::mesh->elements->distribution[t + 1]);
			for (auto i = begin; i != end; ++i) {
				eregions[t][(*i - info::mesh->elements->distribution[t]) * regionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}

	info::mesh->elements->regionMaskSize = regionsBitMaskSize;
	info::mesh->elements->regions = new serializededata<esint, esint>(regionsBitMaskSize, eregions);

	profiler::syncend("fill_region_mask");
	eslog::checkpointln("MESH: REGION MASK FILLED");
}

void computeBoundaryElementsFromNodes(BoundaryRegionStore *bregion, int elementDimension)
{
	profiler::syncstart("compute_boundary_elements_from_nodes");
	if (info::mesh->nodes->elements == NULL) {
		linkNodesAndElements();
	}

	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}

	if (info::mesh->halo->IDs == NULL) {
		exchangeHalo();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<std::pair<esint, esint> > > elements(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = bregion->nodes->begin(t)->begin(); n != bregion->nodes->end(t)->begin(); ++n) {
			auto links = info::mesh->nodes->elements->begin() + *n;
			for (auto e = links->begin(); e != links->end(); ++e) {
				elements[t].push_back(std::make_pair(*e, *n));
			}
		}
	}

	std::vector<size_t> distribution = { 0, elements[0].size() };
	for (size_t t = 1; t < threads; t++) {
		elements[0].insert(elements[0].end(), elements[t].begin(), elements[t].end());
		distribution.push_back(elements[0].size());
	}

	utils::sortWithInplaceMerge(elements[0], distribution);

	esint ebegin = info::mesh->elements->offset;
	esint eend = ebegin + info::mesh->elements->size;

	auto begin = std::lower_bound(elements[0].begin(), elements[0].end(), ebegin,
			[] (const std::pair<esint, esint> &p, esint e) { return p.first < e; });
	auto end = std::lower_bound(elements[0].begin(), elements[0].end(), eend,
			[] (const std::pair<esint, esint> &p, esint e) { return p.first < e; });

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

	int rsize = info::mesh->elements->regionMaskSize;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> nodes, facenodes, lowerElements, lenodes;
		std::vector<esint> tdist, tdata, tcode;
		if (t == 0) {
			tdist.push_back(0);
		}

		int nface;
		esint element, neighbor, prev = 0;
		auto enodes = info::mesh->elements->procNodes->cbegin();
		auto neighbors = info::mesh->elements->faceNeighbors->cbegin();
		const auto &regions = info::mesh->elements->regions->datatarray();

		for (size_t e = tdistribution[t]; e < tdistribution[t + 1]; e++) {
			nodes.push_back((begin + e)->second);
			if ((e + 1 == tdistribution[t + 1] || (begin + e + 1)->first != (begin + e)->first)) {

				element = (begin + e)->first - ebegin;
				utils::sortAndRemoveDuplicates(nodes);

				enodes += element - prev;
				neighbors += element - prev;
				prev = element;

				const auto &fpointers = info::mesh->elements->epointers->datatarray()[element]->facepointers->datatarray();
				auto fnodes = info::mesh->elements->epointers->datatarray()[element]->faces->cbegin();
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
							} else if (element + ebegin < neighbor) {
								if (ebegin <= neighbor && neighbor < eend) {
									neighbor -= ebegin;
									if (memcmp(regions.data() + element * rsize, regions.data() + neighbor * rsize, sizeof(esint) * rsize) != 0) {
										addFace();
									}
								} else {
									neighbor = std::lower_bound(info::mesh->halo->IDs->datatarray().begin(), info::mesh->halo->IDs->datatarray().end(), neighbor) - info::mesh->halo->IDs->datatarray().begin();
									if (memcmp(regions.data() + element * rsize, info::mesh->halo->regions->datatarray().data() + neighbor * rsize, sizeof(esint) * rsize) != 0) {
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

	bregion->procNodes = new serializededata<esint, esint>(edist, edata);
	bregion->epointers = new serializededata<esint, Element*>(1, epointers);
	bregion->dimension = elementDimension;

	profiler::syncend("compute_boundary_elements_from_nodes");
	eslog::checkpoint("MESH: BOUNDARY FROM NODES COMPUTED");
	eslog::param("BOUNDARY", bregion->name.c_str());
	eslog::ln();
}

void synchronizeRegionNodes(std::vector<RegionStore*> &regions)
{
	profiler::syncstart("synchronize_region_nodes");
	std::vector<std::vector<esint> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());
	std::vector<size_t> prevsend(info::mesh->neighbors.size());
	for (auto reg = regions.begin(); reg != regions.end(); ++reg) {
		RegionStore *store = *reg;
		for (size_t n = 0; n < prevsend.size(); ++n) {
			prevsend[n] = sBuffer[n].size();
			sBuffer[n].push_back(0);
		}

		esint prev = 0;
		auto ranks = info::mesh->nodes->ranks->begin();
		for (auto n = store->nodes->datatarray().cbegin(); n != store->nodes->datatarray().cend(); prev = *n++) {
			ranks += *n - prev;

			int rindex = 0;
			for (auto r = ranks->begin(); r != ranks->end(); r++) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbors[rindex] < *r) { ++rindex; }
					sBuffer[rindex].push_back(info::mesh->nodes->IDs->datatarray()[*n]);
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

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::error("ESPRESO internal error: exchange element region nodes.\n");
	}
	profiler::synccheckpoint("exchange");

	std::fill(prevsend.begin(), prevsend.end(), 0);
	std::vector<size_t> prevrecv(info::mesh->neighbors.size());
	for (auto reg = regions.begin(); reg != regions.end(); ++reg) {
		RegionStore *store = *reg;
		std::vector<esint> nnodes;
		for (size_t n = 0; n < rBuffer.size(); ++n) {
			if (info::mesh->neighbors[n] < info::mpi::rank) {
				for (size_t j = prevrecv[n] + 1, s = prevsend[n] + 1; j < prevrecv[n] + rBuffer[n][prevrecv[n]]; ++j) {
					while (s < prevsend[n] + sBuffer[n][prevsend[n]] && sBuffer[n][s] < rBuffer[n][j]) { ++s; }
					if (s == prevsend[n] + sBuffer[n][prevsend[n]] || sBuffer[n][s] != rBuffer[n][j]) {
						auto it = std::find(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().begin() + info::mesh->nodes->uniqInfo.nhalo, rBuffer[n][j]);
						nnodes.push_back(it - info::mesh->nodes->IDs->datatarray().begin());
					}
				}
			} else {
				for (size_t j = prevrecv[n] + 1, s = prevsend[n] + 1; j < prevrecv[n] + rBuffer[n][prevrecv[n]]; ++j) {
					while (s < prevsend[n] + sBuffer[n][prevsend[n]] && sBuffer[n][s] < rBuffer[n][j]) { ++s; }
					if (s == prevsend[n] + sBuffer[n][prevsend[n]] || sBuffer[n][s] != rBuffer[n][j]) {
						auto it = std::lower_bound(info::mesh->nodes->IDs->datatarray().begin() + info::mesh->nodes->uniqInfo.nhalo, info::mesh->nodes->IDs->datatarray().end(), rBuffer[n][j]);
						if (it != info::mesh->nodes->IDs->datatarray().end() && *it == rBuffer[n][j]) {
							nnodes.push_back(it - info::mesh->nodes->IDs->datatarray().begin());
						} else {
							auto it = std::find(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().begin() + info::mesh->nodes->uniqInfo.nhalo, rBuffer[n][j]);
							nnodes.push_back(it - info::mesh->nodes->IDs->datatarray().begin());
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

void computeNodeInfo(std::vector<RegionStore*> &regions)
{
	profiler::syncstart("compute_node_info");
	std::vector<esint> sum, offset;
	for (size_t r = 0; r < regions.size(); ++r) {
		regions[r]->nodeInfo.nhalo = 0;
		for (
				auto n = regions[r]->nodes->datatarray().cbegin();
				n != regions[r]->nodes->datatarray().cend() && *n < info::mesh->nodes->uniqInfo.nhalo;
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

	std::vector<std::vector<esint> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());
	std::vector<size_t> prevsize(info::mesh->neighbors.size()), nsize(info::mesh->neighbors.size());
	std::vector<double> min(3 * regions.size(), std::numeric_limits<double>::max()), max(3 * regions.size(), -std::numeric_limits<double>::max());
	for (size_t r = 0; r < regions.size(); ++r) {
		for (size_t n = 0; n < prevsize.size(); ++n) {
			prevsize[n] = sBuffer[n].size();
			sBuffer[n].push_back(0);
		}

		esint prev = 0, i = 0;
		auto ranks = info::mesh->nodes->ranks->begin();
		for (auto n = regions[r]->nodes->datatarray().cbegin(); n != regions[r]->nodes->datatarray().cend(); prev = *n++, ++i) {
			info::mesh->nodes->coordinates->datatarray()[*n].minmax(min.data() + 3 * r, max.data() + 3 * r);
			ranks += *n - prev;
			if (i < regions[r]->nodeInfo.nhalo) {
				int rindex = 0;
				while (info::mesh->neighbors[rindex] != ranks->front()) { ++rindex; }
				++nsize[rindex];
			} else {
				int rindex = 0;
				for (auto rank = ranks->begin() + 1; rank != ranks->end(); rank++) {
					while (info::mesh->neighbors[rindex] < *rank) { ++rindex; }
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

	if (!Communication::receiveLowerKnownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::error("ESPRESO internal error: receive global offset of a given element region.\n");
	}
	Communication::allReduce(min, Communication::OP::MIN);
	Communication::allReduce(max, Communication::OP::MAX);
	profiler::synccheckpoint("exchange");

	std::fill(prevsize.begin(), prevsize.end(), 0);
	for (size_t r = 0; r < regions.size(); ++r) {
		regions[r]->nodeInfo.min = Point(min[3 * r], min[3 * r + 1], min[3 * r + 2]);
		regions[r]->nodeInfo.max = Point(max[3 * r], max[3 * r + 1], max[3 * r + 2]);
		for (size_t n = 0, begin = 0; n < rBuffer.size(); ++n) {
			if (info::mesh->neighbors[n] < info::mpi::rank) {
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

