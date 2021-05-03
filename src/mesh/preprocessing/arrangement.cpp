
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
}
}

