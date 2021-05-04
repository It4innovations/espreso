
#include "meshpreprocessing.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
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

void partitiate(ElementStore *elements, NodeStore *nodes, ClusterStore *clusters, DomainStore *domains, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, esint parts, bool uniformDecomposition)
{
	profiler::syncstart("partitiate");
	switch (info::ecf->input.decomposition.sequential_decomposer) {
	case DecompositionConfiguration::SequentialDecomposer::NONE:
		break;
	case DecompositionConfiguration::SequentialDecomposer::METIS:
		if (!METIS::islinked() && info::mpi::rank == 0) {
			eslog::warning("ESPRESO run-time event: METIS is not linked, decomposition into domains is skipped.\n");
		}
		break;
	case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
		if (!Scotch::islinked() && info::mpi::rank == 0) {
			eslog::warning("ESPRESO run-time event: SCOTCH is not linked, decomposition into domains is skipped.\n");
		}
		break;
	case DecompositionConfiguration::SequentialDecomposer::KAHIP:
		if (!KaHIP::islinked() && info::mpi::rank == 0) {
			eslog::warning("ESPRESO run-time event: KaHIP is not linked, decomposition into domains is skipped.\n");
		}
		break;
	}

	std::vector<esint> dualDist, dualData;
	computeDecomposedDual(nodes, elements, elementsRegions, neighbors, dualDist, dualData);

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<int> partID(elements->distribution.process.size, -1);

	int nextID = 0;
	for (esint e = 0; e < elements->distribution.process.size; ++e) {
		std::vector<esint> stack;
		if (partID[e] == -1) {
			stack.push_back(e);
			partID[e] = nextID;
			while (stack.size()) {
				esint current = stack.back();
				stack.pop_back();
				for (auto n = dualData.begin() + dualDist[current]; n != dualData.begin() + dualDist[current + 1]; ++n) {
					if (partID[*n] == -1) {
						stack.push_back(*n);
						partID[*n] = nextID;
					}
				}
			}
			nextID++;
		}
	}

	std::vector<int> cluster;
	std::vector<esint> partition(elements->distribution.process.size);

	profiler::synccheckpoint("check_noncontinuity");
	eslog::checkpointln("MESH: CLUSTER NONCONTINUITY CHECKED");

	if (nextID == 1) {
		eslog::checkpointln("MESH: NONCONTINUITY PROCESSED");

		if (uniformDecomposition) {
			esint psize = elements->distribution.process.size / parts;
			for (esint p = 0, offset = 0; p < parts; ++p, offset += psize) {
				std::fill(partition.begin() + offset, partition.begin() + offset + psize, p);
			}
		} else {
			switch (info::ecf->input.decomposition.sequential_decomposer) {
			case DecompositionConfiguration::SequentialDecomposer::NONE:
				break;
			case DecompositionConfiguration::SequentialDecomposer::METIS:
				if (METIS::islinked()) {
					METIS::call(info::ecf->input.decomposition.metis_options,
						elements->distribution.process.size, dualDist.data(), dualData.data(),
						0, NULL, NULL, parts, partition.data());
				}
				profiler::checkpoint("metis");
				break;
			case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
				if (Scotch::islinked()) {
					Scotch::call(info::ecf->input.decomposition.scotch_options,
						elements->distribution.process.size, dualDist.data(), dualData.data(),
						0, NULL, NULL, parts, partition.data());
				}
				profiler::checkpoint("scotch");
				break;
			case DecompositionConfiguration::SequentialDecomposer::KAHIP:
				if (KaHIP::islinked()) {
					KaHIP::call(info::ecf->input.decomposition.kahip_options,
						elements->distribution.process.size, dualDist.data(), dualData.data(),
						0, NULL, NULL, parts, partition.data());
				}
				profiler::checkpoint("kahip");
				break;
			}
		}
		cluster.resize(parts, 0);
		clusters->size = 1;
	} else { // non-continuous dual graph
		// thread x part x elements
		std::vector<std::vector<std::vector<esint> > > tdecomposition(threads, std::vector<std::vector<esint> >(nextID));
		std::vector<std::vector<esint> > tdualsize(threads, std::vector<esint>(nextID));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = elements->distribution.threads[t]; e < elements->distribution.threads[t + 1]; ++e) {
				tdecomposition[t][partID[e]].push_back(e);
				tdualsize[t][partID[e]] += dualDist[e + 1] - dualDist[e];
			}
		}
		std::vector<std::vector<esint> > foffsets(nextID), noffsets(nextID);
		std::vector<esint> partoffset(nextID);
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			foffsets[p].push_back(0);
			noffsets[p].push_back(0);
			for (size_t t = 1; t < threads; t++) {
				foffsets[p].push_back(tdecomposition[0][p].size());
				noffsets[p].push_back(tdualsize[0][p]);
				tdecomposition[0][p].insert(tdecomposition[0][p].end(), tdecomposition[t][p].begin(), tdecomposition[t][p].end());
				tdualsize[0][p] += tdualsize[t][p];
			}
		}
		for (int p = 1; p < nextID; p++) {
			partoffset[p] = partoffset[p - 1] + tdecomposition[0][p - 1].size();
		}

		std::vector<std::vector<esint> > frames(nextID), neighbors(nextID);
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			frames[p].resize(1 + tdecomposition[0][p].size());
			neighbors[p].resize(tdualsize[0][p]);
		}

		// TODO: try parallelization
		// #pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t partindex;
			std::vector<esint> foffset(nextID), noffset(nextID);
			for (int p = 0; p < nextID; p++) {
				foffset[p] = foffsets[p][t];
				noffset[p] = noffsets[p][t];
			}

			for (size_t e = elements->distribution.threads[t]; e < elements->distribution.threads[t + 1]; ++e) {
				partindex = partID[e];

				frames[partindex][++foffset[partindex]] = dualDist[e + 1] - dualDist[e];
				if (e > elements->distribution.threads[t]) {
					frames[partindex][foffset[partindex]] += frames[partindex][foffset[partindex] - 1];
				} else {
					frames[partindex][foffset[partindex]] += noffset[partindex];
				}
				auto node = dualData.begin() + dualDist[e];
				for (esint n = frames[partindex][foffset[partindex]] - (dualDist[e + 1] - dualDist[e]); n < frames[partindex][foffset[partindex]]; ++n, ++node) {
					neighbors[partindex][n] = std::lower_bound(tdecomposition[0][partindex].begin(), tdecomposition[0][partindex].end(), *node) - tdecomposition[0][partindex].begin();
				}
			}
		}

		std::vector<esint> pparts(nextID);

		double averageDomainSize = elements->distribution.process.size / (double)parts;
		size_t partsCounter = 0;
		for (int p = 0; p < nextID; p++) {
			partsCounter += pparts[p] = std::ceil((frames[p].size() - 1) / averageDomainSize);
			cluster.resize(partsCounter, p);
		}
		clusters->size = nextID;

		profiler::checkpoint("process_noncontinuity");
		eslog::checkpointln("MESH: NONCONTINUITY PROCESSED");

		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			switch (info::ecf->input.decomposition.sequential_decomposer) {
			case DecompositionConfiguration::SequentialDecomposer::NONE:
				break;
			case DecompositionConfiguration::SequentialDecomposer::METIS:
				if (METIS::islinked()) {
					METIS::call(info::ecf->input.decomposition.metis_options,
						frames[p].size() - 1, frames[p].data(), neighbors[p].data(),
						0, NULL, NULL, pparts[p], partition.data() + partoffset[p]);
				}
				profiler::checkpoint("metis");
				break;
			case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
				if (Scotch::islinked()) {
					Scotch::call(info::ecf->input.decomposition.scotch_options,
						frames[p].size() - 1, frames[p].data(), neighbors[p].data(),
						0, NULL, NULL, pparts[p], partition.data() + partoffset[p]);
				}
				profiler::checkpoint("scotch");
				break;
			case DecompositionConfiguration::SequentialDecomposer::KAHIP:
				if (KaHIP::islinked()) {
					KaHIP::call(info::ecf->input.decomposition.kahip_options,
						frames[p].size() - 1, frames[p].data(), neighbors[p].data(),
						0, NULL, NULL, pparts[p], partition.data() + partoffset[p]);
				}
				profiler::checkpoint("kahip");
				break;
			}
		}

		std::vector<esint> ppartition = partition;
		nextID = 0;
		for (size_t p = 0; p < tdecomposition[0].size(); p++) {
			for (size_t i = 0; i < tdecomposition[0][p].size(); ++i) {
				partition[tdecomposition[0][p][i]] = ppartition[partoffset[p] + i] + nextID;
			}
			nextID += pparts[p];
		}
	}

	if (uniformDecomposition) {
		eslog::checkpointln("MESH: DOMAINS KEEPED UNIFORM");
	} else {
		switch (info::ecf->input.decomposition.sequential_decomposer) {
		case DecompositionConfiguration::SequentialDecomposer::NONE:
			eslog::checkpointln("MESH: DECOMPOSITION TO DOMAINS SKIPPED"); break;
		case DecompositionConfiguration::SequentialDecomposer::METIS:
			eslog::checkpointln("MESH: DOMAINS COMPUTED BY METIS"); break;
		case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
			eslog::checkpointln("MESH: DOMAINS COMPUTED BY SCOTCH"); break;
		case DecompositionConfiguration::SequentialDecomposer::KAHIP:
			eslog::checkpointln("MESH: DOMAINS COMPUTED BY KAHIP"); break;
		}
	}

	std::vector<esint> permutation(partition.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		if (partition[i] == partition[j]) {
			return i < j;
		}
		return partition[i] < partition[j];
	});

	std::vector<size_t> domainDistribution;
	std::vector<size_t> tdistribution;

	esint partindex = 0;
	auto begin = permutation.begin();
	while (begin != permutation.end()) {
		domainDistribution.push_back(begin - permutation.begin());
		begin = std::lower_bound(begin, permutation.end(), ++partindex, [&] (esint i, esint val) {
			return partition[i] < val;
		});
	}
	domainDistribution.push_back(permutation.size());

	// TODO: improve domain distribution for more complicated decomposition
	if (domainDistribution.size() == threads + 1) {
		tdistribution = std::vector<size_t>(domainDistribution.begin(), domainDistribution.end());
	} else {
		if (domainDistribution.size() < threads + 1) {
			tdistribution = tarray<size_t>::distribute(threads, permutation.size());
		} else {
			double averageThreadSize = elements->distribution.process.size / (double)threads;
			tdistribution.push_back(0);
			for (size_t t = 0; t < threads - 1; t++) {
				auto more = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution.back() + averageThreadSize);
				if (more == domainDistribution.end()) {
					tdistribution.push_back(domainDistribution.back());
				} else {
					auto less = more - 1;
					if (std::fabs(*less - averageThreadSize * (t + 1)) < std::fabs(*more - averageThreadSize * (t + 1))) {
						tdistribution.push_back(*less);
					} else {
						tdistribution.push_back(*more);
					}
				}
			}
			tdistribution.push_back(permutation.size());
		}
	}

	for (size_t i = 1; i < domainDistribution.size(); i++) {
		if (domainDistribution[i - 1] != domainDistribution[i]) {
			domains->cluster.push_back(cluster[i - 1]);
		}
	}
	utils::removeDuplicates(domainDistribution);

	std::vector<size_t> domainCounter(threads);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				++domainCounter[t];
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				++domainCounter[t];
			}
		}
	}

	domains->size = utils::sizesToOffsets(domainCounter);
	domains->offset = domains->size;
	domains->next = domains->size;
	domains->totalSize = Communication::exscan(domains->offset);
	domains->next += domains->offset;
	domainCounter.push_back(domains->size);
	domains->distribution = domainCounter;

	domains->elements.push_back(0);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				domains->elements.push_back(domainDistribution[t + 1]);
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				domains->elements.push_back(*(it + 1));
			}
		}
	}

	profiler::synccheckpoint("arrange_to_domains");
	eslog::checkpointln("MESH: ELEMENTS ARRANGED TO DOMAINS");

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &p = elements->epointers->datatarray();
		for (size_t d = domains->distribution[t]; d < domains->distribution[t + 1]; ++d) {
			std::sort(permutation.begin() + domains->elements[d], permutation.begin() + domains->elements[d + 1], [p] (esint i, esint j) {
				if (p[i]->code != p[j]->code) {
					return p[i]->code < p[j]->code;
				}
				return i < j;
			});
		}
	}
	permuteElements(elements, nodes, elementsRegions, neighbors, permutation, tdistribution);

	profiler::synccheckpoint("arrange_data");
	profiler::syncend("partitiate");
}

void permuteElements(ElementStore *elements, NodeStore *nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors, const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	profiler::syncstart("permute_elements");
	std::vector<esint> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return permutation[i] < permutation[j]; });

	size_t threads = info::env::OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(neighbors.begin(), neighbors.end(), neighbor) - neighbors.begin();
	};

	std::vector<esint> IDBoundaries = Communication::getDistribution(elements->distribution.process.size);
	std::vector<std::vector<std::pair<esint, esint> > > rHalo(neighbors.size());

	if (elements->faceNeighbors != NULL || nodes->elements != NULL) {
		// thread x neighbor x elements(oldID, newID)
		std::vector<std::vector<std::vector<std::pair<esint, esint> > > > sHalo(threads, std::vector<std::vector<std::pair<esint, esint> > >(neighbors.size()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto ranks = nodes->ranks->cbegin(t);
			auto elements = nodes->elements->cbegin(t);
			esint begine = IDBoundaries[info::mpi::rank];
			esint ende   = IDBoundaries[info::mpi::rank + 1];

			for (auto n = nodes->distribution[t]; n < nodes->distribution[t + 1]; ++n, ++ranks, ++elements) {
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != info::mpi::rank) {
						for (auto e = elements->begin(); e != elements->end(); ++e) {
							if (begine <= *e && *e < ende) {
								sHalo[t][n2i(*rank)].push_back(std::make_pair(*e, backpermutation[*e - begine] + begine));
							}
						}
					}
				}
			}

			for (size_t n = 0; n < sHalo[t].size(); ++n) {
				utils::sortAndRemoveDuplicates(sHalo[t][n]);
			}
		}

		utils::mergeThreadedUniqueData(sHalo);

		if (!Communication::exchangeUnknownSize(sHalo[0], rHalo, neighbors)) {
			eslog::internalFailure("exchange halo element new IDs while element permutation.\n");
		}
	}
	profiler::synccheckpoint("prepare");

	auto globalremap = [&] (serializededata<esint, esint>* data, bool sort) {
		if (data == NULL) {
			return;
		}
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			int source;
			for (auto e = data->begin(t); e != data->end(t); ++e) {
				for (auto n = e->begin(); n != e->end(); ++n) {
					if (*n >= 0) {
						source = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), *n + 1) - IDBoundaries.begin() - 1;
						if (source == info::mpi::rank) {
							*n = IDBoundaries[info::mpi::rank] + backpermutation[*n - IDBoundaries[info::mpi::rank]];
						} else {
							*n = std::lower_bound(rHalo[n2i(source)].begin(), rHalo[n2i(source)].end(), std::make_pair(*n, (esint)0))->second;
						}
					}
				}
				if (sort) {
					std::sort(e->begin(), e->end());
				}
			}
		}
	};

	esint firstID = elements->IDs->datatarray().front();
	elements->permute(permutation, distribution);
	std::iota(elements->IDs->datatarray().begin(), elements->IDs->datatarray().end(), firstID);

	globalremap(elements->faceNeighbors, false);
	globalremap(nodes->elements, true);

	for (size_t r = 0; r < elementsRegions.size(); ++r) {
		for (auto n = elementsRegions[r]->elements->datatarray().begin(); n != elementsRegions[r]->elements->datatarray().end(); ++n) {
			*n = backpermutation[*n];
		}
		std::sort(elementsRegions[r]->elements->datatarray().begin(), elementsRegions[r]->elements->datatarray().end());
	}

	profiler::synccheckpoint("remap");
	profiler::syncend("permute_elements");
	eslog::checkpointln("MESH: ELEMENTS PERMUTED");
}

void computeDecomposedDual(NodeStore *nodes, ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors, std::vector<esint> &dualDist, std::vector<esint> &dualData)
{
	profiler::syncstart("compute_decomposed_dual");
	bool separateRegions = info::ecf->input.decomposition.separate_regions;
	bool separateMaterials = info::ecf->input.decomposition.separate_materials;
	bool separateEtypes = info::ecf->input.decomposition.separate_etypes;

	size_t threads = info::env::OMP_NUM_THREADS;
	esint eBegin = elements->distribution.process.offset;
	esint eEnd   = eBegin + elements->distribution.process.size;

	std::vector<esint> dDistribution(elements->distribution.process.size + 1);
	std::vector<std::vector<esint> > dData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;
		int mat1 = 0, mat2 = 0, reg = 0, etype1 = 0, etype2 = 0;
		int rsize = elements->regions->edataSize();

		auto neighs = elements->faceNeighbors->cbegin(t);
		for (size_t e = elements->distribution.threads[t]; e < elements->distribution.threads[t + 1]; ++e, ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1 && eBegin <= *n && *n < eEnd) {
					if (separateMaterials) {
						mat1 = elements->material->datatarray()[e];
						mat2 = elements->material->datatarray()[*n - eBegin];
					}
					if (separateRegions) {
						reg = memcmp(elements->regions->datatarray().data() + e * rsize, elements->regions->datatarray().data() + (*n - eBegin) * rsize, sizeof(esint) * rsize);
					}
					if (separateEtypes) {
						etype1 = (int)elements->epointers->datatarray()[e]->type;
						etype2 = (int)elements->epointers->datatarray()[*n - eBegin]->type;
					}

					if (mat1 == mat2 && !reg && etype1 == etype2) {
						tdata.push_back(*n - eBegin);
					}
				}
			}
			dDistribution[e + 1] = tdata.size();
		}

		dData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dDistribution, elements->distribution.threads);
	for (size_t t = 1; t < threads; t++) {
		dData[0].insert(dData[0].end(), dData[t].begin(), dData[t].end());
	}

	dualDist.swap(dDistribution);
	dualData.swap(dData[0]);

	profiler::syncend("compute_decomposed_dual");
	eslog::checkpointln("MESH: LOCAL DUAL GRAPH COMPUTED");
}

void computeNodeDomainDistribution(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<int> neighborsWithMe)
{
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

	eslog::checkpointln("MESH: NODE TO DOMAINS MAP COMPUTED");
}

void computeLocalIndices(ElementStore *elements, DomainStore *domains)
{
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

	eslog::checkpointln("MESH: LOCAL INDICES COMPUTED");
}

void computeDomainDual(NodeStore *nodes, ElementStore *elements, DomainStore *domains, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe)
{
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

	domains->localDual = new serializededata<esint, esint>(dist, data);
	domains->dual = new serializededata<esint, esint>(distFull, dataFull);

	eslog::checkpointln("MESH: DOMAIN DUAL GRAPH COMPUTED");
}

void computeDomainsSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors)
{
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
				auto faces = epointer->faces->begin();
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

	eslog::checkpointln("MESH: DOMAIN SURFACE COMPUTED");
}

void triangularizeDomainSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors)
{
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

	eslog::checkpointln("MESH: DOMAIN SURFACE TRIANGULARIZED");
}

}
}

