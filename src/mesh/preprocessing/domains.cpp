
#include "meshpreprocessing.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/fetidatastore.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/communication.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "config/ecf/input/decomposition.h"

#include "math/math.h"
#include "wrappers/metis/w.metis.h"

#include "output/visualization/debug.h"

#include <algorithm>
#include <numeric>

namespace espreso {
namespace mesh {

void computeNodeDomainDistribution()
{
	// nID, domain
	std::vector<std::vector<std::pair<esint, esint> > > ntodomains(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<std::pair<esint, esint> > tdata;

		for (esint d = info::mesh->elements->domainDistribution[t]; d != info::mesh->elements->domainDistribution[t + 1]; ++d) {
			std::vector<esint> dnodes(
					(info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d])->begin(),
					(info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d + 1])->begin());

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
	std::vector<std::vector<esint> > rBuffer(info::mesh->neighborsWithMe.size());

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<std::vector<esint> > tBuffer(info::mesh->neighborsWithMe.size());

		auto nranks = info::mesh->nodes->ranks->begin() + info::mesh->nodes->distribution[t];
		auto ntod = std::lower_bound(ntodomains[0].begin(), ntodomains[0].end(), std::pair<esint, esint>(info::mesh->nodes->distribution[t], 0));
		for (size_t n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n, ++nranks) {
			auto begin = ntod;
			while (ntod != ntodomains[0].end() && ntod->first == (esint)n) {
				++ntod;
			}

			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				while (info::mesh->neighborsWithMe[noffset] < *r) {
					++noffset;
				}

				tBuffer[noffset].push_back(ntod - begin);
				for (auto i = begin; i != ntod; ++i) {
					tBuffer[noffset].push_back(info::mesh->elements->firstDomain + i->second);
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

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, info::mesh->neighborsWithMe)) {
		eslog::error("ESPRESO internal error: exchange node domain distribution.\n");
	}

	std::vector<esint> domainsDistrubtion({0});
	std::vector<int> domainsData;

	auto nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> roffset(rBuffer.size());
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		esint noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			while (info::mesh->neighborsWithMe[noffset] < *r) {
				++noffset;
			}

			esint domains = rBuffer[noffset][roffset[noffset]++];
			for (esint d = 0; d < domains; ++d) {
				domainsData.push_back(rBuffer[noffset][roffset[noffset]++]);
			}
		}
		domainsDistrubtion.push_back(domainsData.size());
	}

	std::vector<size_t> ddistribution = info::mesh->nodes->distribution, ddatadistribution = info::mesh->nodes->distribution;
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

	info::mesh->nodes->domains = new serializededata<esint, int>(tarray<esint>(ddistribution, domainsDistrubtion), tarray<int>(ddatadistribution, domainsData));

	eslog::checkpointln("MESH: NODE TO DOMAINS MAP COMPUTED");
}

void computeLocalIndices()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	info::mesh->elements->domainNodes = new serializededata<esint, esint>(*info::mesh->elements->procNodes);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d != info::mesh->elements->domainDistribution[t + 1]; ++d) {
			esint dbegin = (info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d])->begin() - info::mesh->elements->procNodes->datatarray().begin();
			esint dend = (info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d + 1])->begin() - info::mesh->elements->procNodes->datatarray().begin();

			std::vector<esint> dnodes(info::mesh->elements->domainNodes->datatarray().begin() + dbegin, info::mesh->elements->domainNodes->datatarray().begin() + dend);
			utils::sortAndRemoveDuplicates(dnodes);
			for (auto n = info::mesh->elements->domainNodes->datatarray().begin() + dbegin; n != info::mesh->elements->domainNodes->datatarray().begin() + dend; ++n) {
				*n = std::lower_bound(dnodes.begin(), dnodes.end(), *n) - dnodes.begin();
			}
		}
	}

	eslog::checkpointln("MESH: LOCAL INDICES COMPUTED");
}

void computeDomainDual()
{
	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > dist(threads), data(threads);

	esint first = *info::mesh->elements->IDs->datatarray().begin();
	esint last = first + info::mesh->elements->size;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}

		auto neighs = info::mesh->elements->faceNeighbors->cbegin(t);
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			std::vector<esint> ndomains;
			for (esint e = info::mesh->elements->elementsDistribution[d]; e < info::mesh->elements->elementsDistribution[d + 1]; ++e, ++neighs) {
				for (auto n = neighs->begin(); n != neighs->end(); ++n) {
					if (*n != -1 && first <= *n && *n < last) {
						if (*n < first + info::mesh->elements->elementsDistribution[d] || first + info::mesh->elements->elementsDistribution[d + 1] <= *n) {
							esint doffset = 0;
							while (first + info::mesh->elements->elementsDistribution[doffset] <= *n) {
								++doffset;
							}
							ndomains.push_back(doffset - 1);
						}
					}
				}
			}
			utils::sortAndRemoveDuplicates(ndomains);
			tdata.insert(tdata.end(), ndomains.begin(), ndomains.end());
			tdist.push_back(tdata.size());
		}

		dist[t].swap(tdist);
		data[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dist);

	info::mesh->FETIData->domainDual = new serializededata<esint, esint>(dist, data);

	eslog::checkpointln("MESH: DOMAIN DUAL GRAPH COMPUTED");
}

void computeDomainsSurface()
{
	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), ecounter(threads, std::vector<esint>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointer(threads);
	std::vector<std::vector<size_t> > intervals(threads);

	esint eoffset = info::mesh->elements->offset;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tfaces, tfacesDistribution, tecounter((int)Element::CODE::SIZE);
		std::vector<Element*> tfpointer;
		std::vector<size_t> tintervals;
		if (t == 0) {
			tfacesDistribution.push_back(0);
			tintervals.push_back(0);
		}

		auto neighbors = info::mesh->elements->faceNeighbors->cbegin(t);
		auto enodes = info::mesh->elements->procNodes->cbegin(t);
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			esint dbegin = info::mesh->elements->elementsDistribution[d];
			esint dend = info::mesh->elements->elementsDistribution[d + 1];

			for (esint e = dbegin; e < dend; ++e, ++neighbors, ++enodes) {
				auto epointer = info::mesh->elements->epointers->datatarray()[e];
				auto faces = epointer->faces->begin();
				auto facepointer = epointer->facepointers->datatarray().begin();

				for (size_t n = 0; n < neighbors->size(); ++n, ++faces, ++facepointer) {
					if (neighbors->at(n) < dbegin + eoffset || dend + eoffset <= neighbors->at(n)) {
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

	info::mesh->domainsSurface->epointers = new serializededata<esint, Element*>(1, fpointer);
	info::mesh->domainsSurface->ecounters = ecounter[0];

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

	info::mesh->domainsSurface->edistribution = intervals[0];

	std::vector<size_t> tdistribution = { 0 };
	for (size_t t = 0; t < threads; t++) {
		tdistribution.push_back(info::mesh->domainsSurface->edistribution[info::mesh->elements->domainDistribution[t + 1]]);
	}

	if (ecounter[0][(int)Element::CODE::TRIANGLE3] == (esint)info::mesh->domainsSurface->edistribution.back()) {
		serializededata<esint, esint>::balance(3, faces, &tdistribution);
		info::mesh->domainsSurface->enodes = new serializededata<esint, esint>(3, faces);
		info::mesh->domainsSurface->triangles = info::mesh->domainsSurface->enodes;
		info::mesh->domainsSurface->tdistribution = info::mesh->domainsSurface->edistribution;
	} else {
		utils::threadDistributionToFullDistribution(facesDistribution);
		serializededata<esint, esint>::balance(facesDistribution, faces, &tdistribution);
		info::mesh->domainsSurface->enodes = new serializededata<esint, esint>(facesDistribution, faces);
	}

	eslog::checkpointln("MESH: DOMAIN SURFACE COMPUTED");
}

void triangularizeDomainSurface()
{
	if (info::mesh->domainsSurface->enodes == NULL) {
		computeDomainsSurface();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	if (info::mesh->domainsSurface->triangles == NULL) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);

		intervals.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
				auto elements = info::mesh->domainsSurface->enodes->cbegin() + info::mesh->domainsSurface->edistribution[d];
				auto epointer = info::mesh->domainsSurface->epointers->datatarray().begin() + info::mesh->domainsSurface->edistribution[d];

				for (size_t e = info::mesh->domainsSurface->edistribution[d]; e < info::mesh->domainsSurface->edistribution[d + 1]; ++e, ++elements, ++epointer) {
					for (auto n = (*epointer)->triangles->datatarray().cbegin(); n != (*epointer)->triangles->datatarray().cend(); ++n) {
						triangles[t].push_back(elements->at(*n));
					}
				}
				intervals[t].push_back(triangles[t].size() / 3);
			}
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		info::mesh->domainsSurface->tdistribution = intervals[0];
		info::mesh->domainsSurface->triangles = new serializededata<esint, esint>(3, triangles);
	}

	eslog::checkpointln("MESH: DOMAIN SURFACE TRIANGULARIZED");
}

}
}
