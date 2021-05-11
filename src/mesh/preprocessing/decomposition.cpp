
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

static void computeDecomposedDual(const ElementStore *elements, std::vector<esint> &dualDist, std::vector<esint> &dualData)
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

void computeElementsDecomposition(const ElementStore *elements, esint parts, std::vector<size_t> &distribution, std::vector<esint> &permutation)
{
	profiler::syncstart("partitiate");
	bool skip = false;
	switch (info::ecf->input.decomposition.sequential_decomposer) {
	case DecompositionConfiguration::SequentialDecomposer::NONE:
		break;
	case DecompositionConfiguration::SequentialDecomposer::METIS:
		skip = !METIS::islinked();
		if (!METIS::islinked() && info::mpi::rank == 0) {
			eslog::warning("ESPRESO run-time event: METIS is not linked, decomposition into domains is skipped.\n");
		}
		break;
	case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
		skip = !Scotch::islinked();
		if (!Scotch::islinked() && info::mpi::rank == 0) {
			eslog::warning("ESPRESO run-time event: SCOTCH is not linked, decomposition into domains is skipped.\n");
		}
		break;
	case DecompositionConfiguration::SequentialDecomposer::KAHIP:
		skip = !KaHIP::islinked();
		if (!KaHIP::islinked() && info::mpi::rank == 0) {
			eslog::warning("ESPRESO run-time event: KaHIP is not linked, decomposition into domains is skipped.\n");
		}
		break;
	}

	std::vector<esint> dualDist, dualData;
	computeDecomposedDual(elements, dualDist, dualData);

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

	std::vector<esint> partition(partID.begin(), partID.end());

	profiler::synccheckpoint("check_noncontinuity");
	eslog::checkpointln("MESH: CLUSTER NONCONTINUITY CHECKED");

	if (!skip) {
		if (nextID == 1) {
			eslog::checkpointln("MESH: NONCONTINUITY PROCESSED");

			switch (info::ecf->input.decomposition.sequential_decomposer) {
			case DecompositionConfiguration::SequentialDecomposer::NONE:
				break;
			case DecompositionConfiguration::SequentialDecomposer::METIS:
				METIS::call(info::ecf->input.decomposition.metis_options, partition.size(), dualDist.data(), dualData.data(), 0, NULL, NULL, parts, partition.data());
				profiler::checkpoint("metis");
				break;
			case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
				Scotch::call(info::ecf->input.decomposition.scotch_options, partition.size(), dualDist.data(), dualData.data(), 0, NULL, NULL, parts, partition.data());
				profiler::checkpoint("scotch");
				break;
			case DecompositionConfiguration::SequentialDecomposer::KAHIP:
				KaHIP::call(info::ecf->input.decomposition.kahip_options, partition.size(), dualDist.data(), dualData.data(), 0, NULL, NULL, parts, partition.data());
				profiler::checkpoint("kahip");
				break;
			}
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
			profiler::checkpoint("process_noncontinuity");
			eslog::checkpointln("MESH: NONCONTINUITY PROCESSED");

			#pragma omp parallel for
			for (int p = 0; p < nextID; p++) {
				switch (info::ecf->input.decomposition.sequential_decomposer) {
				case DecompositionConfiguration::SequentialDecomposer::NONE:
					break;
				case DecompositionConfiguration::SequentialDecomposer::METIS:
					METIS::call(info::ecf->input.decomposition.metis_options, frames[p].size() - 1, frames[p].data(), neighbors[p].data(), 0, NULL, NULL, pparts[p], partition.data() + partoffset[p]);
					profiler::checkpoint("metis");
					break;
				case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
					Scotch::call(info::ecf->input.decomposition.scotch_options, frames[p].size() - 1, frames[p].data(), neighbors[p].data(), 0, NULL, NULL, pparts[p], partition.data() + partoffset[p]);
					profiler::checkpoint("scotch");
					break;
				case DecompositionConfiguration::SequentialDecomposer::KAHIP:
					KaHIP::call(info::ecf->input.decomposition.kahip_options, frames[p].size() - 1, frames[p].data(), neighbors[p].data(), 0, NULL, NULL, pparts[p], partition.data() + partoffset[p]);
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

	permutation.clear();
	permutation.resize(partition.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		if (partition[i] == partition[j]) {
			if (elements->epointers->datatarray()[i]->code == elements->epointers->datatarray()[j]->code) {
				return i < j;
			}
			return elements->epointers->datatarray()[i]->code < elements->epointers->datatarray()[j]->code;
		}
		return partition[i] < partition[j];
	});

	distribution.clear();
	esint partindex = 0;
	auto begin = permutation.begin();
	while (begin != permutation.end()) {
		distribution.push_back(begin - permutation.begin());
		begin = std::lower_bound(begin, permutation.end(), ++partindex, [&] (esint i, esint val) { return partition[i] < val; });
	}
	distribution.push_back(permutation.size());
	profiler::syncend("partitiate");
	eslog::checkpointln("MESH: ELEMENT PERMUTATION COMPUTED");
}

void permuteElements(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces, std::vector<int> &neighbors, std::vector<size_t> &distribution, const std::vector<esint> &permutation)
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

	std::vector<size_t> tdistribution;
	// TODO: improve domain distribution for more complicated decomposition
	if (distribution.size() == threads + 1) {
		tdistribution = std::vector<size_t>(distribution.begin(), distribution.end());
	} else {
		if (distribution.size() < threads + 1) {
			tdistribution = tarray<size_t>::distribute(threads, permutation.size());
		} else {
			double averageThreadSize = elements->distribution.process.size / (double)threads;
			tdistribution.push_back(0);
			for (size_t t = 0; t < threads - 1; t++) {
				auto more = std::lower_bound(distribution.begin(), distribution.end(), tdistribution.back() + averageThreadSize);
				if (more == distribution.end()) {
					tdistribution.push_back(distribution.back());
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

	utils::removeDuplicates(distribution);

	std::vector<size_t> domainCounter(threads);
	for (size_t t = 0; t < threads; t++) {
		if (distribution.size() < threads + 1) {
			if (t < distribution.size() - 1) {
				++domainCounter[t];
			}
		} else {
			auto begin = std::lower_bound(distribution.begin(), distribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(distribution.begin(), distribution.end(), tdistribution[t + 1]);
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

	domains->elements.clear();
	domains->elements.push_back(0);
	for (size_t t = 0; t < threads; t++) {
		if (distribution.size() < threads + 1) {
			if (t < distribution.size() - 1) {
				domains->elements.push_back(distribution[t + 1]);
			}
		} else {
			auto begin = std::lower_bound(distribution.begin(), distribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(distribution.begin(), distribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				domains->elements.push_back(*(it + 1));
			}
		}
	}

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
	elements->permute(permutation, tdistribution);
	std::iota(elements->IDs->datatarray().begin(), elements->IDs->datatarray().end(), firstID);

	globalremap(elements->faceNeighbors, false);
	globalremap(nodes->elements, true);

	for (size_t r = 0; r < elementsRegions.size(); ++r) {
		for (auto n = elementsRegions[r]->elements->datatarray().begin(); n != elementsRegions[r]->elements->datatarray().end(); ++n) {
			*n = backpermutation[*n];
		}
		std::sort(elementsRegions[r]->elements->datatarray().begin(), elementsRegions[r]->elements->datatarray().end());
	}

	for (size_t r = 0; r < boundaryRegions.size(); ++r) {
		if (boundaryRegions[r]->emembership != NULL) {
			for (auto em = boundaryRegions[r]->emembership->begin(); em != boundaryRegions[r]->emembership->end(); ++em) {
				em->at(0) = backpermutation[em->at(0)];
			}
		}
	}

	for (size_t r = 0; r < contactInterfaces.size(); ++r) {
		if (contactInterfaces[r]->emembership != NULL) {
			for (auto em = contactInterfaces[r]->emembership->begin(); em != contactInterfaces[r]->emembership->end(); ++em) {
				em->at(0) = backpermutation[em->at(0)];
			}
		}
	}

	profiler::synccheckpoint("remap");
	profiler::syncend("permute_elements");
	eslog::checkpointln("MESH: ELEMENTS PERMUTED");
}

} // namespace mesh
} // namespace espreso
