
#include "meshpreprocessing.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/contactinterfacestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/fetidatastore.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "config/ecf/input/decomposition.h"

#include "math/math.h"
#include "wrappers/metis/w.metis.h"
#include "wrappers/scotch/w.scotch.h"
#include "wrappers/kahip/w.kahip.h"

#include "output/visualization/debug.h"

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
				if (memcmp(elements->regions->datatarray().data() + e * elements->regions->edataSize(), elements->regions->datatarray().data() + (e - 1) * elements->regions->edataSize(), sizeof(esint) * elements->regions->edataSize())) {
					iboundaries[t].push_back(e);
				}
			}
		}
	}
	utils::mergeThreadedUniqueData(iboundaries);
	profiler::synccheckpoint("iboundaries");

	auto addregions = [&] (esint i) {
		auto regs = (elements->regions->begin() + i)->data();
		std::vector<int> regions;
		for (size_t r = 1; r < elements->regions->edataSize() * sizeof(esint); ++r) {
			esint maskOffset = r / (8 * sizeof(esint));
			esint bit = 1 << (r % (8 * sizeof(esint)));
			if (regs[maskOffset] & bit) {
				regions.push_back(r);
			}
		}
		if (regions.size() == 0) {
			elements->eintervals.back().region = 0;
		}
		if (regions.size() == 1) {
			elements->eintervals.back().region = regions.front();
		}
		if (regions.size() > 1) {
			elements->eintervals.back().region = -1;
			elements->eintervals.back().regions = regions;
		}
	};

	elements->eintervals.clear();
	elements->eintervals.push_back(ElementsInterval(0, 0));
	elements->eintervals.back().domain = domains->offset;
	elements->eintervals.back().code = static_cast<int>(elements->epointers->datatarray().front()->code);
	addregions(0);
	elements->eintervalsDistribution.push_back(0);
	for (size_t i = 0; i < iboundaries[0].size(); i++) {
		elements->eintervals.back().end = iboundaries[0][i];
		elements->eintervals.push_back(ElementsInterval(iboundaries[0][i], iboundaries[0][i]));
		const std::vector<esint> &edist = domains->elements;
		elements->eintervals.back().domain = std::lower_bound(edist.begin(), edist.end(), elements->eintervals.back().begin + 1) - edist.begin() - 1 + domains->offset;
		elements->eintervals.back().code = static_cast<int>(elements->epointers->datatarray()[elements->eintervals.back().begin]->code);
		addregions(elements->eintervals.back().begin);
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

void computeRegionsBoundaryIntervals(const ElementStore *elements, const DomainStore *domains, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces)
{
	profiler::syncstart("arrange_boudary_regions");
	int threads = info::env::OMP_NUM_THREADS;

	std::vector<BoundaryRegionStore*> allRegions;
	allRegions.insert(allRegions.end(), boundaryRegions.begin(), boundaryRegions.end());
	allRegions.insert(allRegions.end(), contactInterfaces.begin(), contactInterfaces.end());

	for (size_t r = 0; r < allRegions.size(); r++) {
		BoundaryRegionStore *store = allRegions[r];
		if (store->dimension) {
			std::vector<esint> eregion(store->distribution.process.size);
			#pragma omp parallel for
			for (int t = 0; t < threads; t++) {
				auto region = eregion.begin() + store->emembership->datatarray().distribution()[t] / 2;
				for (auto parent = store->emembership->begin(t); parent != store->emembership->end(t); ++parent, ++region) {
					*region = std::lower_bound(elements->eintervals.begin(), elements->eintervals.end(), parent->at(0) + 1, [&] (const ElementsInterval &ei, esint e) { return ei.end < e; }) - elements->eintervals.begin();
				}
			}

			std::vector<esint> permutation(store->elements->structures());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
				if (eregion[i] == eregion[j]) {
					if (store->epointers->datatarray()[i]->code == store->epointers->datatarray()[j]->code) {
						return store->emembership->datatarray()[2 * i] < store->emembership->datatarray()[2 * j];
					}
					return store->epointers->datatarray()[i]->code < store->epointers->datatarray()[j]->code;
				}
				return eregion[i] < eregion[j];
			});

			store->eintervals.clear();
			if (store->distribution.process.size) {
				std::vector<size_t> distribution = { 0 };
				auto dend = domains->distribution.begin() + 1;
				esint begin = 0;
				for (auto prev = permutation.begin(), currect = prev + 1; currect != permutation.end(); prev = currect++) {
					if (eregion[*prev] != eregion[*currect] || store->epointers->datatarray()[*prev]->code != store->epointers->datatarray()[*currect]->code) {
						store->eintervals.push_back(ElementsInterval(begin, currect - permutation.begin()));
						store->eintervals.back().code = static_cast<int>(store->epointers->datatarray()[*prev]->code);
						store->eintervals.back().domain = elements->eintervals[eregion[*prev]].domain;
						begin = currect - permutation.begin();
					}
					if (eregion[*prev] != eregion[*currect]) {
						while (elements->eintervals[eregion[*currect]].domain >= (esint)*dend) {
							distribution.push_back(currect - permutation.begin());
							++dend;
						}
					}
				}
				distribution.resize(threads + 1, eregion.size());
				store->eintervals.push_back(ElementsInterval(begin, eregion.size()));
				store->eintervals.back().code = static_cast<int>(store->epointers->datatarray()[permutation.back()]->code);
				store->eintervals.back().domain = elements->eintervals[eregion[permutation.back()]].domain;
				store->permute(permutation, distribution);
			}
			store->eintervalsDistribution.resize(domains->size + 1);
			auto ei = store->eintervals.begin();
			for (esint d = 0; d < domains->size; ++d) {
				while (ei != store->eintervals.end() && ei->domain <= d) {
					++ei;
				}
				store->eintervalsDistribution[d + 1] = ei - store->eintervals.begin();
			}
		}
	}

	profiler::syncend("arrange_boudary_regions");
	eslog::checkpointln("MESH: BOUNDARY REGIONS ARRANGED");
}

void computeNodeDomainDistribution(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<int> neighborsWithMe)
{
	profiler::syncstart("compute_node_domain_distribution");
	// nID, domain
	std::vector<std::vector<std::pair<esint, esint> > > ntodomains(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<std::pair<esint, esint> > tdata;

		for (size_t d = domains->distribution[t]; d != domains->distribution[t + 1]; ++d) {
			std::vector<esint> dnodes(
					(elements->nodes->begin() + domains->elements[d])->begin(),
					(elements->nodes->begin() + domains->elements[d + 1])->begin());

			utils::sortAndRemoveDuplicates(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				tdata.push_back(std::pair<esint, esint>(dnodes[i], d));
			}
		}

		ntodomains[t].swap(tdata);
	}

	for (int t = 1; t < info::env::OMP_NUM_THREADS; t++) {
		ntodomains[0].insert(ntodomains[0].end(), ntodomains[t].begin(), ntodomains[t].end());
	}

	std::sort(ntodomains[0].begin(), ntodomains[0].end());

	std::vector<std::vector<std::vector<esint> > > sBuffer(info::env::OMP_NUM_THREADS);
	std::vector<std::vector<esint> > rBuffer(neighborsWithMe.size());

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<std::vector<esint> > tBuffer(neighborsWithMe.size());

		auto nranks = nodes->ranks->begin() + nodes->distribution[t];
		auto ntod = std::lower_bound(ntodomains[0].begin(), ntodomains[0].end(), std::pair<esint, esint>(nodes->distribution[t], 0));
		for (size_t n = nodes->distribution[t]; n < nodes->distribution[t + 1]; ++n, ++nranks) {
			auto begin = ntod;
			while (ntod != ntodomains[0].end() && ntod->first == (esint)n) {
				++ntod;
			}

			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				while (neighborsWithMe[noffset] < *r) {
					++noffset;
				}

				tBuffer[noffset].push_back(ntod - begin);
				for (auto i = begin; i != ntod; ++i) {
					tBuffer[noffset].push_back(domains->offset + i->second);
				}
			}
		}
		sBuffer[t].swap(tBuffer);
	}

	for (int t = 1; t < info::env::OMP_NUM_THREADS; t++) {
		for (size_t n = 0; n < sBuffer[t].size(); n++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, neighborsWithMe)) {
		eslog::internalFailure("exchange node domain distribution.\n");
	}

	std::vector<esint> domainsDistrubtion({0});
	std::vector<int> domainsData;

	auto nranks = nodes->ranks->begin();
	std::vector<esint> roffset(rBuffer.size());
	for (esint n = 0; n < nodes->size; ++n, ++nranks) {
		esint noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			while (neighborsWithMe[noffset] < *r) {
				++noffset;
			}

			esint domains = rBuffer[noffset][roffset[noffset]++];
			for (esint d = 0; d < domains; ++d) {
				domainsData.push_back(rBuffer[noffset][roffset[noffset]++]);
			}
		}
		domainsDistrubtion.push_back(domainsData.size());
	}

	std::vector<size_t> ddistribution = nodes->distribution, ddatadistribution = nodes->distribution;
	for (int t = 1; t < info::env::OMP_NUM_THREADS; t++) {
		++ddistribution[t];
		if (ddistribution[t] < domainsDistrubtion.size()) {
			ddatadistribution[t] = domainsDistrubtion[ddistribution[t]];
		} else {
			ddatadistribution[t] = domainsDistrubtion[ddistribution[info::env::OMP_NUM_THREADS] - 1];
		}
	}
	++ddistribution[info::env::OMP_NUM_THREADS];
	ddatadistribution[info::env::OMP_NUM_THREADS] = domainsDistrubtion[ddistribution[info::env::OMP_NUM_THREADS] - 1];

	nodes->domains = new serializededata<esint, int>(tarray<esint>(ddistribution, domainsDistrubtion), tarray<int>(ddatadistribution, domainsData));

	profiler::syncend("compute_node_domain_distribution");
	eslog::checkpointln("MESH: NODE TO DOMAINS MAP COMPUTED");
}

void computeLocalIndices(ElementStore *elements, DomainStore *domains)
{
	profiler::syncstart("compute_local_indices");
	size_t threads = info::env::OMP_NUM_THREADS;

	domains->nodes = new serializededata<esint, esint>(*elements->nodes);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = domains->distribution[t]; d != domains->distribution[t + 1]; ++d) {
			esint dbegin = (elements->nodes->begin() + domains->elements[d])->begin() - elements->nodes->datatarray().begin();
			esint dend = (elements->nodes->begin() + domains->elements[d + 1])->begin() - elements->nodes->datatarray().begin();

			std::vector<esint> dnodes(domains->nodes->datatarray().begin() + dbegin, domains->nodes->datatarray().begin() + dend);
			utils::sortAndRemoveDuplicates(dnodes);
			for (auto n = domains->nodes->datatarray().begin() + dbegin; n != domains->nodes->datatarray().begin() + dend; ++n) {
				*n = std::lower_bound(dnodes.begin(), dnodes.end(), *n) - dnodes.begin();
			}
		}
	}

	profiler::syncend("compute_local_indices");
	eslog::checkpointln("MESH: LOCAL INDICES COMPUTED");
}

void computeDomainDual(NodeStore *nodes, ElementStore *elements, DomainStore *domains, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe)
{
	profiler::syncstart("compute_domain_dual");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::pair<esint, esint> > sBuffer, gBuffer;
	std::vector<std::vector<std::pair<esint, esint> > > rBuffer(neighborsWithMe.size());

	for (esint d = 0; d < domains->size; ++d) {
		sBuffer.push_back(std::make_pair(elements->distribution.process.offset + domains->elements[d + 1], d + domains->offset));
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighborsWithMe)) {
		eslog::internalFailure("cannot exchange distribution of elements to domains.\n");
	}

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		gBuffer.insert(gBuffer.end(), rBuffer[n].begin(), rBuffer[n].end());
	}

	std::vector<std::vector<esint> > dist(threads), data(threads);
	std::vector<std::vector<esint> > distFull(threads), dataFull(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdist, tdata, tdistFull, tdataFull;
		if (t == 0) {
			tdist.push_back(0);
			tdistFull.push_back(0);
		}

		auto neighs = elements->faceNeighbors->cbegin(t);
		for (size_t d = domains->distribution[t]; d < domains->distribution[t + 1]; d++) {
			std::vector<esint> ndomains, ndomainsFull;
			for (esint e = domains->elements[d]; e < domains->elements[d + 1]; ++e, ++neighs) {
				for (auto n = neighs->begin(); n != neighs->end(); ++n) {
					if (*n != -1) {
						if (*n < elements->distribution.process.offset + domains->elements[d] || elements->distribution.process.offset + domains->elements[d + 1] <= *n) {
							auto it = std::lower_bound(gBuffer.begin(), gBuffer.end(), *n, [] (const std::pair<esint, esint> &info, const esint &e) { return info.first <= e; });
							ndomainsFull.push_back(it->second);
							if (domains->isLocal(ndomainsFull.back())) {
								ndomains.push_back(ndomainsFull.back() - domains->offset);
							}
						}
					}
				}
			}
			utils::sortAndRemoveDuplicates(ndomains);
			utils::sortAndRemoveDuplicates(ndomainsFull);
			tdata.insert(tdata.end(), ndomains.begin(), ndomains.end());
			tdist.push_back(tdata.size());
			tdataFull.insert(tdataFull.end(), ndomainsFull.begin(), ndomainsFull.end());
			tdistFull.push_back(tdataFull.size());
		}

		dist[t].swap(tdist);
		data[t].swap(tdata);
		distFull[t].swap(tdistFull);
		dataFull[t].swap(tdataFull);
	}

	utils::threadDistributionToFullDistribution(dist);
	utils::threadDistributionToFullDistribution(distFull);

	if (domains->localDual) delete domains->localDual;
	if (domains->dual) delete domains->dual;
	domains->localDual = new serializededata<esint, esint>(dist, data);
	domains->dual = new serializededata<esint, esint>(distFull, dataFull);

	profiler::syncend("compute_domain_dual");
	eslog::checkpointln("MESH: DOMAIN DUAL GRAPH COMPUTED");
}

void computeClustersDistribution(DomainStore *domains, ClusterStore *clusters)
{
	profiler::syncstart("compute_cluster_distribution");
	domains->cluster.clear();
	domains->cluster.resize(domains->size, -1);

	int cluster = 0;
	for (esint d = 0; d < domains->size; ++d) {
		std::vector<esint> stack;
		if (domains->cluster[d] == -1) {
			stack.push_back(d);
			domains->cluster[d] = cluster;
			while (stack.size()) {
				auto neighs = domains->localDual->begin() + stack.back();
				stack.pop_back();
				for (auto n = neighs->begin(); n != neighs->end(); ++n) {
					if (domains->cluster[*n] == -1) {
						stack.push_back(*n);
						domains->cluster[*n] = cluster;
					}
				}
			}
			cluster++;
		}
	}
	clusters->offset = clusters->next = clusters->size = cluster;
	clusters->totalSize = Communication::exscan(clusters->offset);
	clusters->next += clusters->size;
	profiler::syncend("compute_cluster_distribution");
	eslog::checkpointln("MESH: CLUSTERS DISTRIBUTION COMPUTED");
}

void computeDomainsSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors)
{
	profiler::syncend("compute_domains_surface");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), ecounter(threads, std::vector<esint>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointer(threads);
	std::vector<std::vector<size_t> > intervals(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tfaces, tfacesDistribution, tecounter((int)Element::CODE::SIZE);
		std::vector<Element*> tfpointer;
		std::vector<size_t> tintervals;
		if (t == 0) {
			tfacesDistribution.push_back(0);
			tintervals.push_back(0);
		}

		auto neighbors = elements->faceNeighbors->cbegin(t);
		auto enodes = elements->nodes->cbegin(t);
		for (size_t d = domains->distribution[t]; d < domains->distribution[t + 1]; d++) {
			esint dbegin = domains->elements[d];
			esint dend = domains->elements[d + 1];

			for (esint e = dbegin; e < dend; ++e, ++neighbors, ++enodes) {
				auto epointer = elements->epointers->datatarray()[e];
				auto faces = epointer->faceList->begin();
				auto facepointer = epointer->facepointers->datatarray().begin();

				for (size_t n = 0; n < neighbors->size(); ++n, ++faces, ++facepointer) {
					if (neighbors->at(n) < dbegin + elements->distribution.process.offset || dend + elements->distribution.process.offset <= neighbors->at(n)) {
						for (auto f = faces->begin(); f != faces->end(); ++f) {
							tfaces.push_back(enodes->at(*f));
						}
						tfacesDistribution.push_back(tfaces.size());
						tfpointer.push_back(*facepointer);
						++tecounter[(int)(*facepointer)->code];
					}
				}
			}
			tintervals.push_back(tfacesDistribution.size());
		}

		faces[t].swap(tfaces);
		facesDistribution[t].swap(tfacesDistribution);
		fpointer[t].swap(tfpointer);
		ecounter[t].swap(tecounter);
		intervals[t].swap(tintervals);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t e = 0; e < ecounter[0].size(); e++) {
			ecounter[0][e] += ecounter[t][e];
		}
	}

	domainsSurface->epointers = new serializededata<esint, Element*>(1, fpointer);
	domainsSurface->ecounters = ecounter[0];

	for (size_t i = 1; i < intervals[0].size(); i++) {
		--intervals[0][i];
	}
	esint tsize = facesDistribution[0].size() - 1;
	for (size_t t = 1; t < threads; t++) {
		for (size_t i = 0; i < intervals[t].size(); i++) {
			intervals[t][i] += tsize;
		}
		tsize += facesDistribution[t].size();
	}
	utils::mergeThreadedUniqueData(intervals);
	utils::sortAndRemoveDuplicates(intervals[0]);

	domainsSurface->edistribution = intervals[0];

	std::vector<size_t> tdistribution = { 0 };
	for (size_t t = 0; t < threads; t++) {
		tdistribution.push_back(domainsSurface->edistribution[domains->distribution[t + 1]]);
	}

	if (ecounter[0][(int)Element::CODE::TRIANGLE3] == (esint)domainsSurface->edistribution.back()) {
		serializededata<esint, esint>::balance(3, faces, &tdistribution);
		domainsSurface->enodes = new serializededata<esint, esint>(3, faces);
		domainsSurface->triangles = domainsSurface->enodes;
		domainsSurface->tdistribution = domainsSurface->edistribution;
	} else {
		utils::threadDistributionToFullDistribution(facesDistribution);
		serializededata<esint, esint>::balance(facesDistribution, faces, &tdistribution);
		domainsSurface->enodes = new serializededata<esint, esint>(facesDistribution, faces);
	}

	profiler::syncend("compute_domains_surface");
	eslog::checkpointln("MESH: DOMAIN SURFACE COMPUTED");
}

void triangularizeDomainSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors)
{
	profiler::syncstart("triangulate_domains_surface");
	if (domainsSurface->enodes == NULL) {
		computeDomainsSurface(nodes, elements, domains, domainsSurface, neighbors);
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	if (domainsSurface->triangles == NULL) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);

		intervals.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t d = domains->distribution[t]; d < domains->distribution[t + 1]; d++) {
				auto elements = domainsSurface->enodes->cbegin() + domainsSurface->edistribution[d];
				auto epointer = domainsSurface->epointers->datatarray().begin() + domainsSurface->edistribution[d];

				for (size_t e = domainsSurface->edistribution[d]; e < domainsSurface->edistribution[d + 1]; ++e, ++elements, ++epointer) {
					for (auto n = (*epointer)->triangles->datatarray().cbegin(); n != (*epointer)->triangles->datatarray().cend(); ++n) {
						triangles[t].push_back(elements->at(*n));
					}
				}
				intervals[t].push_back(triangles[t].size() / 3);
			}
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		domainsSurface->tdistribution = intervals[0];
		domainsSurface->triangles = new serializededata<esint, esint>(3, triangles);
	}

	profiler::syncend("triangulate_domains_surface");
	eslog::checkpointln("MESH: DOMAIN SURFACE TRIANGULARIZED");
}

}
}

