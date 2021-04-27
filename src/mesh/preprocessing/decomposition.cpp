
#include "meshpreprocessing.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"

#include "config/ecf/input/decomposition.h"

#include "output/visualization/debug.h"

#include "wrappers/kahip/w.kahip.h"
#include "wrappers/metis/w.metis.h"
#include "wrappers/parmetis/w.parmetis.h"
#include "wrappers/ptscotch/w.ptscotch.h"
#include "wrappers/scotch/w.scotch.h"

#include <algorithm>
#include <numeric>

namespace espreso {
namespace mesh {

esint getSFCDecomposition(std::vector<esint> &partition)
{
	profiler::syncstart("get_sfc_decomposition");
	int threads = info::env::OMP_NUM_THREADS;

	if (info::mesh->elements->centers == NULL) {
		computeElementsCenters();
	}

	HilbertCurve sfc(info::mesh->dimension, SFCDEPTH, info::mesh->nodes->coordinates->datatarray().size(), info::mesh->nodes->coordinates->datatarray().data());

	std::vector<esint> buckets(info::mesh->elements->process.size), borders;
	std::vector<esint> permutation(info::mesh->elements->process.size);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		auto center = info::mesh->elements->centers->datatarray().cbegin(t);
		for (size_t e = info::mesh->elements->threading[t]; e != info::mesh->elements->threading[t + 1]; ++e, ++center) {
			buckets[e] = sfc.getBucket(*center);
		}
	}

	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		return buckets[i] < buckets[j];
	});
	profiler::synccheckpoint("sfc_buckets");

	if (!Communication::computeSplitters(buckets, permutation, borders)) {
		eslog::internalFailure("cannot compute splitters.\n");
	}
	borders.back() = sfc.buckets(sfc.depth());
//	Communication::computeSFCBalancedBorders(sfc, buckets, permutation, borders);
	profiler::synccheckpoint("compute_splitters");

	if (info::mesh->elements->process.size) {
		auto border = std::lower_bound(borders.begin(), borders.end(), buckets[permutation[0]] + 1);
		for (esint e = 0; e != info::mesh->elements->process.size; ++e) {
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

	for (esint e = 0; e < info::mesh->elements->process.size; ++e) {
		dcenters[partition[e]] += info::mesh->elements->centers->datatarray()[e];
		++dsize[partition[e]];
	}

	Communication::allReduce(dcenters.data(), sumcenters.data(), 3 * info::mpi::size, MPI_DOUBLE, MPI_SUM);
	Communication::allReduce(dsize.data(), sumsize.data(), info::mpi::size, MPITools::getType<esint>().mpitype, MPI_SUM);

	for (int r = 0; r < info::mpi::size; ++r) {
		dcenters[r] = sumcenters[r] / sumsize[r];
	}

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		for (size_t e = info::mesh->elements->threading[t]; e < info::mesh->elements->threading[t + 1]; ++e) {
			Point &center = info::mesh->elements->centers->datatarray()[e];
			for (int r = 0; r < info::mpi::size; ++r) {
				if ((dcenters[r] - center).length() < (dcenters[partition[e]] - center).length()) {
					partition[e] = r;
				}
			}
		}
	}

	return 0; // edge cut is not computed
}

esint callParallelDecomposer(std::vector<esint> &eframes, std::vector<esint> &eneighbors, std::vector<esint> &partition)
{
	DebugOutput::meshDual(eframes, eneighbors);

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
		edgecut = getSFCDecomposition(partition);
		eslog::checkpointln("MESH: RECLUSTERIZED ACCORDING SFC");
		return edgecut;
	}

	profiler::syncstart("call_parallel_decomposer");
	if (info::mpi::size <= info::ecf->input.third_party_scalability_limit && info::ecf->input.decomposition.parallel_decomposer != DecompositionConfiguration::ParallelDecomposer::METIS) {
		std::vector<esint> edistribution = Communication::getDistribution<esint>(partition.size(), &MPITools::subset->across);

		switch (info::ecf->input.decomposition.parallel_decomposer) {
		case DecompositionConfiguration::ParallelDecomposer::NONE: break;
		case DecompositionConfiguration::ParallelDecomposer::METIS: break; // never accessed
		case DecompositionConfiguration::ParallelDecomposer::PARMETIS:
			edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
					MPITools::subset->across, edistribution.data(), eframes.data(), eneighbors.data(),
					0, NULL, 0, NULL, NULL, partition.data());

			profiler::synccheckpoint("parmetis");
			eslog::checkpointln("MESH: RECLUSTERIZED BY PARMETIS");

			if (info::ecf->input.decomposition.parmetis_options.refinement) {
				esint prev = 2 * edgecut;
				while (1.01 * edgecut < prev) {
					prev = edgecut;
					edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_RefineKway,
							MPITools::subset->across, edistribution.data(), eframes.data(), eneighbors.data(),
							0, NULL, 0, NULL, NULL, partition.data());
				}
				profiler::synccheckpoint("parmetis");
				eslog::checkpointln("MESH: CLUSTERIZATION REFINED BY PARMETIS");
			}
			break;
		case DecompositionConfiguration::ParallelDecomposer::PTSCOTCH:
			edgecut = PTScotch::call(
					MPITools::subset->across, edistribution.data(), eframes.data(), eneighbors.data(),
					0, NULL, 0, NULL, NULL, partition.data());

			profiler::synccheckpoint("ptscotch");
			eslog::checkpointln("MESH: RECLUSTERIZED BY PTSCOTCH");
			break;
		case DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE: break; // never accessed
		}
	} else {
		MPIType type = MPITools::getType<esint>();
		std::vector<esint> gframes, gneighbors, gpartition, edistribution;
		std::vector<size_t> offsets;

		if (info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::METIS) {
			Communication::gatherUnknownSize(eframes, gframes, &MPITools::singleton->within);
			Communication::gatherUnknownSize(eneighbors, gneighbors, &MPITools::singleton->within);
			Communication::gatherUnknownSize(partition, gpartition, offsets, &MPITools::singleton->within);
		} else {
			Communication::gatherUnknownSize(eframes, gframes, &MPITools::subset->within);
			Communication::gatherUnknownSize(eneighbors, gneighbors, &MPITools::subset->within);
			Communication::gatherUnknownSize(partition, gpartition, offsets, &MPITools::subset->within);
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
			if (MPITools::subset->within.rank == 0) {
				fixframes();
				edistribution = Communication::getDistribution<esint>(gpartition.size(), &MPITools::subset->across);

				switch (info::ecf->input.decomposition.parallel_decomposer) {
				case DecompositionConfiguration::ParallelDecomposer::NONE: break;
				case DecompositionConfiguration::ParallelDecomposer::METIS: break;// never accessed
				case DecompositionConfiguration::ParallelDecomposer::PARMETIS:
					edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
							MPITools::subset->across, edistribution.data(), gframes.data(), gneighbors.data(),
							0, NULL, 0, NULL, NULL, gpartition.data());
					profiler::checkpoint("parmetis");
					break;
				case DecompositionConfiguration::ParallelDecomposer::PTSCOTCH:
					edgecut = PTScotch::call(
							MPITools::subset->across, edistribution.data(), gframes.data(), gneighbors.data(),
							0, NULL, 0, NULL, NULL, gpartition.data());
					profiler::checkpoint("ptscotch");
					break;
				case DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE: break; // never accessed
				}
			}

			profiler::synccheckpoint("synchronize");
			Communication::barrier(&MPITools::subset->within);
		}

		switch (info::ecf->input.decomposition.parallel_decomposer) {
		case DecompositionConfiguration::ParallelDecomposer::NONE: eslog::checkpointln("MESH: RECLUSTERIZED SKIPPED"); break;
		case DecompositionConfiguration::ParallelDecomposer::METIS: eslog::checkpointln("MESH: RECLUSTERIZED BY METIS"); break;
		case DecompositionConfiguration::ParallelDecomposer::PARMETIS: eslog::checkpointln("MESH: RECLUSTERIZED BY PARMETIS"); break;
		case DecompositionConfiguration::ParallelDecomposer::PTSCOTCH: eslog::checkpointln("MESH: RECLUSTERIZED BY PTSCOTCH"); break;
		case DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE: break; // never accessed
		}

		if (
				info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::PARMETIS &&
				info::ecf->input.decomposition.parmetis_options.refinement) {

			if (MPITools::subset->within.rank == 0) {
				esint prev = 2 * edgecut;
				while (1.01 * edgecut < prev) {
					prev = edgecut;
					edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_RefineKway,
							MPITools::subset->across, edistribution.data(), gframes.data(), gneighbors.data(),
							0, NULL, 0, NULL, NULL, gpartition.data());
				}
			}
			Communication::barrier(&MPITools::subset->within);
			profiler::synccheckpoint("refinement");
			eslog::checkpointln("MESH: CLUSTERS REFINED BY PARMETIS");
		}

		if (info::ecf->input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::METIS) {
			Communication::broadcastUnknownSize(offsets, &MPITools::singleton->within);
			Communication::scatterv(gpartition, partition, offsets, &MPITools::singleton->within);
			Communication::broadcast(&edgecut, 1, type.mpitype, 0, &MPITools::singleton->within);
		} else {
			Communication::broadcastUnknownSize(offsets, &MPITools::subset->within);
			Communication::scatterv(gpartition, partition, offsets, &MPITools::subset->within);
			Communication::broadcast(&edgecut, 1, type.mpitype, 0, &MPITools::subset->within);
		}


		profiler::synccheckpoint("expand");
		eslog::checkpointln("MESH: MPI PROCESSES EXPANDED");
	}

	profiler::syncend("call_parallel_decomposer");
	return edgecut;
}

void reclusterize()
{
	if (info::mpi::size == 1) {
		return;
	}

	profiler::syncstart("reclusterize");
	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}

	if (info::mesh->halo->IDs == NULL) {
		exchangeHalo();
	}

	// Disable due to horible scalability
//	if (info::mesh->elements->centers == NULL) {
//		computeElementsCenters();
//	}

	bool separateRegions = info::ecf->input.decomposition.separate_regions;
	bool separateMaterials = info::ecf->input.decomposition.separate_materials;
	bool separateEtypes = info::ecf->input.decomposition.separate_etypes;

	if (separateRegions && info::mesh->elements->regions == NULL) {
		fillRegionMask();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	esint eoffset = info::mesh->elements->process.offset;

	std::vector<esint> dDistribution(info::mesh->elements->process.size + 1);
	std::vector<std::vector<esint> > dData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;
		int mat1 = 0, mat2 = 0, reg = 0, etype1 = 0, etype2 = 0;
		int rsize = info::mesh->elements->regionMaskSize;
		esint hindex = 0;

		auto neighs = info::mesh->elements->faceNeighbors->cbegin(t);
		for (size_t e = info::mesh->elements->threading[t]; e < info::mesh->elements->threading[t + 1]; ++e, ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1) {
					if ((separateRegions || separateMaterials || separateEtypes) && (*n < eoffset || eoffset + info::mesh->elements->process.size <= *n)) {
						hindex = std::lower_bound(info::mesh->halo->IDs->datatarray().begin(), info::mesh->halo->IDs->datatarray().end(), *n) - info::mesh->halo->IDs->datatarray().begin();
					}
					if (separateMaterials) {
						mat1 = info::mesh->elements->material->datatarray()[e];
						if (*n < eoffset || eoffset + info::mesh->elements->process.size <= *n) {
							mat2 = info::mesh->halo->material->datatarray()[hindex];
						} else {
							mat2 = info::mesh->elements->material->datatarray()[*n - eoffset];
						}
					}
					if (separateRegions) {
						if (*n < eoffset || eoffset + info::mesh->elements->process.size <= *n) {
							reg = memcmp(info::mesh->elements->regions->datatarray().data() + e * rsize, info::mesh->halo->regions->datatarray().data() + hindex * rsize, sizeof(esint) * rsize);
						} else {
							reg = memcmp(info::mesh->elements->regions->datatarray().data() + e * rsize, info::mesh->elements->regions->datatarray().data() + (*n - eoffset) * rsize, sizeof(esint) * rsize);
						}
					}
					if (separateEtypes) {
						etype1 = (int)info::mesh->elements->epointers->datatarray()[e]->type;
						if (*n < eoffset || eoffset + info::mesh->elements->process.size <= *n) {
							etype2 = (int)info::mesh->halo->epointers->datatarray()[hindex]->type;
						} else {
							etype2 = (int)info::mesh->elements->epointers->datatarray()[*n - eoffset]->type;
						}
					}

					if (mat1 == mat2 && !reg && etype1 == etype2) {
						tdata.push_back(*n);
					}
				}
			}
			dDistribution[e + 1] = tdata.size();
		}

		dData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dDistribution, info::mesh->elements->threading);
	for (size_t t = 1; t < threads; t++) {
		dData[0].insert(dData[0].end(), dData[t].begin(), dData[t].end());
	}

	std::vector<esint> partition(info::mesh->elements->process.size, info::mpi::rank);

	profiler::synccheckpoint("compute_dual");
	eslog::checkpointln("MESH: DUAL GRAPH COMPUTED");

	callParallelDecomposer(dDistribution, dData.front(), partition);
	profiler::synccheckpoint("decompose");

	if (info::ecf->input.decomposition.parallel_decomposer != DecompositionConfiguration::ParallelDecomposer::NONE) {
		exchangeElements(partition);
	}
	profiler::synccheckpoint("exchange");
	profiler::syncend("reclusterize");
}

void partitiate(esint parts, bool uniformDecomposition)
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
	computeDecomposedDual(dualDist, dualData);

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<int> partID(info::mesh->elements->process.size, -1);

	int nextID = 0;
	for (esint e = 0; e < info::mesh->elements->process.size; ++e) {
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

	std::vector<int> clusters;
	std::vector<esint> partition(info::mesh->elements->process.size);

	profiler::synccheckpoint("check_noncontinuity");
	eslog::checkpointln("MESH: CLUSTER NONCONTINUITY CHECKED");

	if (nextID == 1) {
		eslog::checkpointln("MESH: NONCONTINUITY PROCESSED");

		if (uniformDecomposition) {
			esint psize = info::mesh->elements->process.size / parts;
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
						info::mesh->elements->process.size, dualDist.data(), dualData.data(),
						0, NULL, NULL, parts, partition.data());
				}
				profiler::checkpoint("metis");
				break;
			case DecompositionConfiguration::SequentialDecomposer::SCOTCH:
				if (Scotch::islinked()) {
					Scotch::call(info::ecf->input.decomposition.scotch_options,
						info::mesh->elements->process.size, dualDist.data(), dualData.data(),
						0, NULL, NULL, parts, partition.data());
				}
				profiler::checkpoint("scotch");
				break;
			case DecompositionConfiguration::SequentialDecomposer::KAHIP:
				if (KaHIP::islinked()) {
					KaHIP::call(info::ecf->input.decomposition.kahip_options,
						info::mesh->elements->process.size, dualDist.data(), dualData.data(),
						0, NULL, NULL, parts, partition.data());
				}
				profiler::checkpoint("kahip");
				break;
			}
		}
		clusters.resize(parts, 0);
		info::mesh->clusters->size = 1;
	} else { // non-continuous dual graph
		// thread x part x elements
		std::vector<std::vector<std::vector<esint> > > tdecomposition(threads, std::vector<std::vector<esint> >(nextID));
		std::vector<std::vector<esint> > tdualsize(threads, std::vector<esint>(nextID));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = info::mesh->elements->threading[t]; e < info::mesh->elements->threading[t + 1]; ++e) {
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

			for (size_t e = info::mesh->elements->threading[t]; e < info::mesh->elements->threading[t + 1]; ++e) {
				partindex = partID[e];

				frames[partindex][++foffset[partindex]] = dualDist[e + 1] - dualDist[e];
				if (e > info::mesh->elements->threading[t]) {
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

		double averageDomainSize = info::mesh->elements->process.size / (double)parts;
		size_t partsCounter = 0;
		for (int p = 0; p < nextID; p++) {
			partsCounter += pparts[p] = std::ceil((frames[p].size() - 1) / averageDomainSize);
			clusters.resize(partsCounter, p);
		}
		info::mesh->clusters->size = nextID;

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
			double averageThreadSize = info::mesh->elements->process.size / (double)threads;
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
			info::mesh->domains->cluster.push_back(clusters[i - 1]);
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

	info::mesh->domains->size = utils::sizesToOffsets(domainCounter);
	info::mesh->domains->offset = info::mesh->domains->size;
	Communication::exscan(info::mesh->domains->offset);
	domainCounter.push_back(info::mesh->domains->size);
	info::mesh->domains->distribution = domainCounter;

	info::mesh->domains->elements.push_back(0);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				info::mesh->domains->elements.push_back(domainDistribution[t + 1]);
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				info::mesh->domains->elements.push_back(*(it + 1));
			}
		}
	}

	profiler::synccheckpoint("arrange_to_domains");
	eslog::checkpointln("MESH: ELEMENTS ARRANGED TO DOMAINS");

	arrangeElementsPermutation(permutation);
	permuteElements(permutation, tdistribution);

	arrangeNodes();
	profiler::synccheckpoint("arrange_data");
	profiler::syncend("partitiate");
}

void exchangeElements(const std::vector<esint> &partition)
{
	profiler::syncstart("exchange_elements");
	if (info::mesh->nodes->elements == NULL) {
		// need for correctly update nodes ranks
		linkNodesAndElements();
	}

	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}

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
			for (size_t e = info::mesh->elements->threading[t]; e < info::mesh->elements->threading[t + 1]; ++e) {
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

	ElementStore *elements = new ElementStore();

	std::vector<std::vector<esint> >  elemsIDs(threads);
	std::vector<std::vector<int> >      elemsBody(threads);
	std::vector<std::vector<int> >      elemsMaterial(threads);
	std::vector<std::vector<Element*> > elemsEpointer(threads);
	std::vector<std::vector<esint> >  elemsNodesDistribution(threads);
	std::vector<std::vector<esint> >  elemsNodesData(threads);
	std::vector<std::vector<esint> >  elemsNeighborsDistribution(threads);
	std::vector<std::vector<esint> >  elemsNeighborsData(threads);
	std::vector<std::vector<esint> >  elemsRegions(threads);

	NodeStore *nodes = new NodeStore();

	std::vector<std::vector<esint> >  nodesIDs(threads);
	std::vector<std::vector<Point> >    nodesCoordinates(threads);
	std::vector<std::vector<esint> >  nodesElemsDistribution(threads);
	std::vector<std::vector<esint> >  nodesElemsData(threads);
	std::vector<std::vector<esint> >  nodesRegions(threads);

	std::vector<std::vector<std::vector<esint> > >  boundaryEDistribution(info::mesh->boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<esint> > >  boundaryEData(info::mesh->boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<Element*> > > boundaryEPointers(info::mesh->boundaryRegions.size(), std::vector<std::vector<Element*> >(threads));

	// regions are transfered via mask
	int eregionsBitMaskSize = info::mesh->elementsRegions.size() / (8 * sizeof(esint)) + (info::mesh->elementsRegions.size() % (8 * sizeof(esint)) ? 1 : 0);
	int bregionsBitMaskSize = info::mesh->boundaryRegions.size() / (8 * sizeof(esint)) + (info::mesh->boundaryRegions.size() % (8 * sizeof(esint)) ? 1 : 0);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, body, material, code, dualsize, dualdata, nodesize, nodeindices)
	std::vector<std::vector<std::vector<esint> > > sElements(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rElements;

	// threads x target x nodes(id, point, linksize, links, regionMask) + size
	std::vector<std::vector<std::vector<esint> > > sNodes(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rNodes;

	// threads x target x boundary(prefix, (code, nodes))
	std::vector<std::vector<std::vector<esint> > > sBoundary(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>(info::mesh->boundaryRegions.size())));
	std::vector<std::vector<esint> > rBoundary;

	// Step 1: Serialize element data

	std::vector<esint> regionElementMask(info::mesh->elements->process.size * eregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));
			auto begin = std::lower_bound(info::mesh->elementsRegions[r]->elements->datatarray().begin(), info::mesh->elementsRegions[r]->elements->datatarray().end(), info::mesh->elements->threading[t]);
			auto end = std::lower_bound(info::mesh->elementsRegions[r]->elements->datatarray().begin(), info::mesh->elementsRegions[r]->elements->datatarray().end(), info::mesh->elements->threading[t + 1]);
			for (auto i = begin; i != end; ++i) {
				regionElementMask[*i * eregionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto IDs = info::mesh->elements->IDs->datatarray().data();
		auto body = info::mesh->elements->body->datatarray().data();
		auto material = info::mesh->elements->material->datatarray().data();
		auto epointer = info::mesh->elements->epointers->datatarray().data();
		auto enodes = info::mesh->elements->nodes->cbegin(t);
		auto eneighbors = info::mesh->elements->faceNeighbors->cbegin(t);
		auto nIDs = info::mesh->nodes->IDs->datatarray().data();

		std::vector<std::vector<esint> > tsElements(targets.size(), std::vector<esint>({ 0 }));

		std::vector<esint>  telemsIDs;
		std::vector<int>      telemsBody;
		std::vector<int>      telemsMaterial;
		std::vector<Element*> telemsEpointer;
		std::vector<esint>  telemsNodesDistribution;
		std::vector<esint>  telemsNodesData;
		std::vector<esint>  telemsNeighborsDistribution;
		std::vector<esint>  telemsNeighborsData;
		std::vector<esint>  telemsRegions;
		if (t == 0) {
			telemsNodesDistribution.push_back(0);
			telemsNeighborsDistribution.push_back(0);
		}

		// estimation
		telemsIDs.reserve(1.5 * info::mesh->elements->process.size / threads);
		telemsBody.reserve(1.5 * info::mesh->elements->process.size / threads);
		telemsMaterial.reserve(1.5 * info::mesh->elements->process.size / threads);
		telemsEpointer.reserve(1.5 * info::mesh->elements->process.size / threads);
		telemsNodesDistribution.reserve(1.5 * info::mesh->elements->process.size / threads);
		telemsNeighborsDistribution.reserve(1.5 * info::mesh->elements->process.size / threads);
		telemsRegions.reserve(1.5 * info::mesh->elements->process.size / threads);

		size_t target;
		for (size_t e = info::mesh->elements->threading[t]; e < info::mesh->elements->threading[t + 1]; ++e, ++enodes, ++eneighbors) {
			if (partition[e] == info::mpi::rank) {
				telemsIDs.push_back(IDs[e]);
				telemsBody.push_back(body[e]);
				telemsMaterial.push_back(material[e]);
				telemsEpointer.push_back(epointer[e]);
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					telemsNodesData.push_back(nIDs[*n]);
				}
				telemsNodesDistribution.push_back(telemsNodesData.size());
				telemsNeighborsData.insert(telemsNeighborsData.end(), eneighbors->begin(), eneighbors->end());
				telemsNeighborsDistribution.push_back(telemsNeighborsData.size());
				telemsRegions.insert(telemsRegions.end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			} else {
				target = t2i(partition[e]);
				tsElements[target].insert(tsElements[target].end(), { IDs[e], body[e], material[e], static_cast<int>(epointer[e]->code) });
				tsElements[target].push_back(enodes->size());
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					tsElements[target].push_back(nIDs[*n]);
				}
				tsElements[target].push_back(eneighbors->size());
				tsElements[target].insert(tsElements[target].end(), eneighbors->begin(), eneighbors->end());
				tsElements[target].insert(tsElements[target].end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			}
		}

		elemsIDs[t].swap(telemsIDs);
		elemsBody[t].swap(telemsBody);
		elemsMaterial[t].swap(telemsMaterial);
		elemsEpointer[t].swap(telemsEpointer);
		elemsNodesDistribution[t].swap(telemsNodesDistribution);
		elemsNodesData[t].swap(telemsNodesData);
		elemsNeighborsDistribution[t].swap(telemsNeighborsDistribution);
		elemsNeighborsData[t].swap(telemsNeighborsData);
		elemsRegions[t].swap(telemsRegions);

		sElements[t].swap(tsElements);
	}

	profiler::synccheckpoint("serialize_elements");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE ELEMENTS");

	// Step 2: Serialize node data

	std::vector<esint> regionNodeMask(info::mesh->nodes->size * bregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->nodes) {
				maskOffset = r / (8 * sizeof(esint));
				esint bit = (esint)1 << (r % (8 * sizeof(esint)));
				auto begin = std::lower_bound(info::mesh->boundaryRegions[r]->nodes->datatarray().begin(), info::mesh->boundaryRegions[r]->nodes->datatarray().end(), info::mesh->nodes->distribution[t]);
				auto end = std::lower_bound(info::mesh->boundaryRegions[r]->nodes->datatarray().begin(), info::mesh->boundaryRegions[r]->nodes->datatarray().end(), info::mesh->nodes->distribution[t + 1]);
				for (auto i = begin; i != end; ++i) {
					regionNodeMask[*i * bregionsBitMaskSize + maskOffset] |= bit;
				}
			}
		}
	}

	esint eBegin = Communication::getDistribution(info::mesh->elements->process.size)[info::mpi::rank];
	esint eEnd = eBegin + info::mesh->elements->process.size;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = info::mesh->nodes->IDs->datatarray();
		const auto &coordinates = info::mesh->nodes->coordinates->datatarray();
		auto elems = info::mesh->nodes->elements->cbegin(t);

		std::vector<esint>  tnodesIDs;
		std::vector<Point>    tnodesCoordinates;
		std::vector<esint>  tnodesElemsDistribution;
		std::vector<esint>  tnodesElemsData;
		std::vector<esint>  tnodesRegions;

		if (t == 0) {
			tnodesElemsDistribution.push_back(0);
		}

		std::vector<std::vector<esint> > tsNodes(targets.size(), std::vector<esint>({ 0 }));

		tnodesIDs.reserve(1.5 * info::mesh->nodes->size / threads);
		tnodesCoordinates.reserve(1.5 * info::mesh->nodes->size / threads);
		tnodesElemsDistribution.reserve(1.5 * info::mesh->nodes->size / threads);
		tnodesRegions.reserve(1.5 * info::mesh->nodes->size / threads);

		size_t target;
		std::vector<bool> last(targets.size() + 1); // targets + me
		for (size_t n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n, ++elems) {
			std::fill(last.begin(), last.end(), false);
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				if (eBegin <= *e && *e < eEnd) {
					target = t2i(partition[*e - eBegin]);
					if (!last[target] && partition[*e - eBegin] != info::mpi::rank) {
						tsNodes[target].push_back(IDs[n]);
						tsNodes[target].insert(tsNodes[target].end(), reinterpret_cast<const esint*>(coordinates.data() + n), reinterpret_cast<const esint*>(coordinates.data() + n + 1));
						tsNodes[target].push_back(elems->size());
						tsNodes[target].insert(tsNodes[target].end(), elems->begin(), elems->end());
						tsNodes[target].insert(tsNodes[target].end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last[target] = true;
					}
					if (!last.back() && partition[*e - eBegin] == info::mpi::rank) {
						tnodesIDs.push_back(IDs[n]);
						tnodesCoordinates.push_back(coordinates[n]);
						tnodesElemsData.insert(tnodesElemsData.end(), elems->begin(), elems->end());
						tnodesElemsDistribution.push_back(tnodesElemsData.size());
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
		nodesRegions[t].swap(tnodesRegions);

		sNodes[t].swap(tsNodes);
	}

	profiler::synccheckpoint("serialize_nodes");
	eslog::checkpointln("EXCHANGE EL: SERIALIZE NODES");

	// Step 2.1: Serialize boundary regions data

	std::vector<esint> emembership;

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			emembership.clear();
			emembership.resize(info::mesh->boundaryRegions[r]->distribution.back());
			std::vector<size_t> distribution = info::mesh->boundaryRegions[r]->distribution;

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + distribution[t];
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
							if (counter == enodes->size() && eBegin <= nlinks[i]) {
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
				auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + distribution[t];
				const auto &IDs = info::mesh->nodes->IDs->datatarray();
				const auto &epointer = info::mesh->boundaryRegions[r]->epointers->datatarray();
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					if (partition[emembership[e]] == info::mpi::rank) {
						tboundaryEPointers.push_back(epointer[e]);
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							tboundaryEData.push_back(IDs[*n]);
						}
						tboundaryEDistribution.push_back(tboundaryEData.size());
					} else {
						target = t2i(partition[emembership[e]]);
						sBoundary[t][target].push_back(static_cast<int>(epointer[e]->code));
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							sBoundary[t][target].push_back(IDs[*n]);
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

			std::vector<esint>  telemsIDs;
			std::vector<int>      telemsBody;
			std::vector<int>      telemsMaterial;
			std::vector<Element*> telemsEpointer;
			std::vector<esint>  telemsNodesDistribution;
			std::vector<esint>  telemsNodesData;
			std::vector<esint>  telemsNeighborsDistribution;
			std::vector<esint>  telemsNeighborsData;
			std::vector<esint>  telemsRegions;

			telemsIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsBody.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsMaterial.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsEpointer.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsNodesDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsNeighborsDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsRegions.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0, neighOffset = 0;

			if (elemsNodesDistribution[t].size()) {
				distOffset = elemsNodesDistribution[t].back();
			}
			if (elemsNeighborsDistribution[t].size()) {
				neighOffset = elemsNeighborsDistribution[t].back();
			}
			if (t == 0 && elemsNodesDistribution[t].size() == 0) {
				telemsNodesDistribution.push_back(0);
			}
			if (t == 0 && elemsNeighborsDistribution[t].size() == 0) {
				telemsNeighborsDistribution.push_back(0);
			}

			for (size_t e = rdistribution[t] + 1; e < rdistribution[t + 1]; ) {
				telemsIDs.push_back(rElements[i][e++]);
				telemsBody.push_back(rElements[i][e++]);
				telemsMaterial.push_back(rElements[i][e++]);
				telemsEpointer.push_back(&Mesh::edata[rElements[i][e++]]);
				telemsNodesData.insert(telemsNodesData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNodesDistribution.push_back(telemsNodesData.size() + distOffset);
				e += rElements[i][e++]; // nodes + nodes size
				telemsNeighborsData.insert(telemsNeighborsData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNeighborsDistribution.push_back(telemsNeighborsData.size() + neighOffset);
				e += rElements[i][e++]; // neighbors + neighbors size
				telemsRegions.insert(telemsRegions.end(), rElements[i].begin() + e, rElements[i].begin() + e + eregionsBitMaskSize);
				e += eregionsBitMaskSize;
			}

			elemsIDs[t].insert(elemsIDs[t].end(), telemsIDs.begin(), telemsIDs.end());
			elemsBody[t].insert(elemsBody[t].end(), telemsBody.begin(), telemsBody.end());
			elemsMaterial[t].insert(elemsMaterial[t].end(), telemsMaterial.begin(), telemsMaterial.end());
			elemsEpointer[t].insert(elemsEpointer[t].end(), telemsEpointer.begin(), telemsEpointer.end());
			elemsNodesDistribution[t].insert(elemsNodesDistribution[t].end(), telemsNodesDistribution.begin(), telemsNodesDistribution.end());
			elemsNodesData[t].insert(elemsNodesData[t].end(), telemsNodesData.begin(), telemsNodesData.end());
			elemsNeighborsDistribution[t].insert(elemsNeighborsDistribution[t].end(), telemsNeighborsDistribution.begin(), telemsNeighborsDistribution.end());
			elemsNeighborsData[t].insert(elemsNeighborsData[t].end(), telemsNeighborsData.begin(), telemsNeighborsData.end());
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
			std::vector<esint>  tnodesRegions;
			std::vector<esint>  tnodeSet;
			std::vector<esint>  tnpermutation;

			tnodesIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesCoordinates.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesElemsDistribution.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesRegions.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodeSet.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0;

			if (nodesElemsDistribution[t].size()) {
				distOffset = nodesElemsDistribution[t].back();
			}
			if (t == 0 && nodesElemsDistribution[t].size() == 0) {
				tnodesElemsDistribution.push_back(0);
			}

			auto nodesetit = nodeset.begin();
			if (rdistribution[t] + 1 < rdistribution[t + 1]) {
				nodesetit = std::lower_bound(nodeset.begin(), nodeset.end(), rNodes[i][rdistribution[t] + 1]);
			}
			for (size_t n = rdistribution[t] + 1; n < rdistribution[t + 1]; ) {
				tnpermutation.push_back(n);
				n += 1 + sizeof(Point) / sizeof(esint); // id, Point
				n += 1 + rNodes[i][n]; // linksize, links
				n += bregionsBitMaskSize; // region mask
			}
			std::sort(tnpermutation.begin(), tnpermutation.end(), [&] (esint n1, esint n2) {
				return rNodes[i][n1] < rNodes[i][n2];
			});
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
					tnodesRegions.insert(tnodesRegions.end(), rNodes[i].begin() + index, rNodes[i].begin() + index + bregionsBitMaskSize);
					index += bregionsBitMaskSize; // region mask
				}
			}

			nodesIDs[t].insert(nodesIDs[t].end(), tnodesIDs.begin(), tnodesIDs.end());
			nodesCoordinates[t].insert(nodesCoordinates[t].end(), tnodesCoordinates.begin(), tnodesCoordinates.end());
			nodesElemsDistribution[t].insert(nodesElemsDistribution[t].end(), tnodesElemsDistribution.begin(), tnodesElemsDistribution.end());
			nodesElemsData[t].insert(nodesElemsData[t].end(), tnodesElemsData.begin(), tnodesElemsData.end());
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
	utils::threadDistributionToFullDistribution(nodesElemsDistribution);

	profiler::synccheckpoint("deserialize_nodes");
	eslog::checkpointln("EXCHANGE EL: DESERIALIZE NODES");

	// Step 4: Deserialize boundary data
	for (size_t n = 0; n < rBoundary.size(); ++n) {
		std::vector<std::vector<esint> > toffset(info::mesh->boundaryRegions.size()), tsize(info::mesh->boundaryRegions.size());
		esint offset = 0, p = 0;
		for (size_t t = 0; t < threads; t++) {
			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
				toffset[r].push_back(offset + info::mesh->boundaryRegions.size());
				tsize[r].push_back(rBoundary[n][p + r]);
				offset += rBoundary[n][p + r];
			}
			offset += info::mesh->boundaryRegions.size();
			p = offset;
		}
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
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
					tboundaryEData.insert(tboundaryEData.end(), rBoundary[n].begin() + i, rBoundary[n].begin() + i + tboundaryEPointers.back()->nodes);
					tboundaryEDistribution.push_back(tboundaryEData.size() + distOffset);
					i += tboundaryEPointers.back()->nodes;
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
		elemDistribution[t] = elemDistribution[t - 1] + elemsIDs[t - 1].size();
	}

	elements->IDs = new serializededata<esint, esint>(1, elemsIDs);
	elements->body = new serializededata<esint, int>(1, elemsBody);
	elements->material = new serializededata<esint, int>(1, elemsMaterial);
	elements->epointers = new serializededata<esint, Element*>(1, elemsEpointer);
	elements->nodes = new serializededata<esint, esint>(elemsNodesDistribution, elemsNodesData); // global IDs

	elements->regionMaskSize = eregionsBitMaskSize;
	elements->regions = new serializededata<esint, esint>(eregionsBitMaskSize, elemsRegions);

	elements->faceNeighbors = new serializededata<esint, esint>(elemsNeighborsDistribution, elemsNeighborsData);

	elements->process.size = elements->IDs->structures();
	elements->threading = elements->IDs->datatarray().distribution();

	// Step 5: Balance node data to threads
	std::vector<size_t> nodeDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		nodeDistribution[t] = nodeDistribution[t - 1] + nodesIDs[t - 1].size();
	}

	serializededata<esint, esint>::balance(1, nodesIDs);
	serializededata<esint, Point>::balance(1, nodesCoordinates);
	serializededata<esint, esint>::balance(nodesElemsDistribution, nodesElemsData);

	nodes->IDs = new serializededata<esint, esint>(1, nodesIDs);
	nodes->coordinates = new serializededata<esint, Point>(1, nodesCoordinates);
	nodes->elements = new serializededata<esint, esint>(nodesElemsDistribution, nodesElemsData);
	nodes->size = nodes->IDs->datatarray().size();
	nodes->distribution = nodes->IDs->datatarray().distribution();

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			delete info::mesh->boundaryRegions[r]->procNodes;
			delete info::mesh->boundaryRegions[r]->epointers;

			utils::threadDistributionToFullDistribution(boundaryEDistribution[r]);
			info::mesh->boundaryRegions[r]->procNodes = new serializededata<esint, esint>(boundaryEDistribution[r], boundaryEData[r]);
			info::mesh->boundaryRegions[r]->epointers = new serializededata<esint, Element*>(1, boundaryEPointers[r]);
		}
	}

	for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
		esint maskOffset = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));
		delete info::mesh->elementsRegions[r]->elements;
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
		info::mesh->elementsRegions[r]->elements = new serializededata<esint, esint>(1, regionelems);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->nodes) {
			esint maskOffset = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));
			delete info::mesh->boundaryRegions[r]->nodes;
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
			info::mesh->boundaryRegions[r]->nodes = new serializededata<esint, esint>(1, regionnodes);
		}
	}

	std::vector<esint> eIDsOLD = Communication::getDistribution(info::mesh->elements->process.size);
	std::vector<esint> eIDsNEW = Communication::getDistribution(elements->process.size);

	for (size_t t = 1; t < threads; ++t) {
		elemsIDs[0].insert(elemsIDs[0].end(), elemsIDs[t].begin(), elemsIDs[t].end());
	}

	std::vector<esint> epermutation(elemsIDs[0].size());
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return elemsIDs[0][i] < elemsIDs[0][j]; });

	std::vector<esint> sortedElements(epermutation.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		for (size_t e = elemDistribution[t]; e < elemDistribution[t + 1]; e++) {
			sortedElements[e] = elemsIDs[0][epermutation[e]];
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
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
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
		for (auto neighbors = elements->faceNeighbors->begin(t); neighbors != elements->faceNeighbors->end(t); ++neighbors) {
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


	std::vector<esint> permutation(nodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return nodes->IDs->datatarray()[i] < nodes->IDs->datatarray()[j]; });

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = elements->nodes->begin(t)->begin(); n != elements->nodes->end(t)->begin(); ++n) {
			*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (esint i, esint val) {
				return nodes->IDs->datatarray()[i] < val;
			});
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->procNodes) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = info::mesh->boundaryRegions[r]->procNodes->begin(t)->begin(); n != info::mesh->boundaryRegions[r]->procNodes->end(t)->begin(); ++n) {
					*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (esint i, esint val) {
						return nodes->IDs->datatarray()[i] < val;
					});
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
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
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

	nodes->ranks = new serializededata<esint, int>(rankBoundaries, rankData);

	std::iota(elements->IDs->datatarray().begin(), elements->IDs->datatarray().end(), eIDsNEW[info::mpi::rank]);
	elements->process.offset = eIDsNEW[info::mpi::rank];
	elements->process.totalSize = info::mesh->elements->process.totalSize;
	std::swap(info::mesh->elements, elements);
	std::swap(info::mesh->nodes, nodes);
	delete info::mesh->halo;
	info::mesh->halo = new ElementStore();
	info::mesh->neighbors.clear();
	for (size_t t = 0; t < IDtargets[0].size(); t++) {
		if (IDtargets[0][t] != info::mpi::rank) {
			info::mesh->neighbors.push_back(IDtargets[0][t]);
		}
	}
	info::mesh->neighborsWithMe = info::mesh->neighbors;
	info::mesh->neighborsWithMe.push_back(info::mpi::rank);
	std::sort(info::mesh->neighborsWithMe.begin(), info::mesh->neighborsWithMe.end());

	delete elements;
	delete nodes;

	profiler::synccheckpoint("finish");
	profiler::syncend("exchange_elements");
	eslog::endln("EXCHANGE EL: FINISH");
	eslog::checkpointln("MESH: ELEMENTS EXCHANGED");
}

void permuteElements(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	profiler::syncstart("permute_elements");
	if (info::mesh->nodes->elements == NULL) {
		linkNodesAndElements();
	}

	std::vector<esint> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return permutation[i] < permutation[j]; });

	size_t threads = info::env::OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(info::mesh->neighbors.begin(), info::mesh->neighbors.end(), neighbor) - info::mesh->neighbors.begin();
	};

	std::vector<esint> IDBoundaries = Communication::getDistribution(info::mesh->elements->process.size);
	std::vector<std::vector<std::pair<esint, esint> > > rHalo(info::mesh->neighbors.size());

	if (info::mesh->elements->faceNeighbors != NULL || info::mesh->nodes->elements != NULL) {
		// thread x neighbor x elements(oldID, newID)
		std::vector<std::vector<std::vector<std::pair<esint, esint> > > > sHalo(threads, std::vector<std::vector<std::pair<esint, esint> > >(info::mesh->neighbors.size()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto ranks = info::mesh->nodes->ranks->cbegin(t);
			auto elements = info::mesh->nodes->elements->cbegin(t);
			esint begine = IDBoundaries[info::mpi::rank];
			esint ende   = IDBoundaries[info::mpi::rank + 1];

			for (auto n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n, ++ranks, ++elements) {
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

		if (!Communication::exchangeUnknownSize(sHalo[0], rHalo, info::mesh->neighbors)) {
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

	esint firstID = info::mesh->elements->IDs->datatarray().front();
	info::mesh->elements->permute(permutation, distribution);
	std::iota(info::mesh->elements->IDs->datatarray().begin(), info::mesh->elements->IDs->datatarray().end(), firstID);

	globalremap(info::mesh->elements->faceNeighbors, false);
	globalremap(info::mesh->nodes->elements, true);

	for (size_t r = 0; r < info::mesh->elementsRegions.size(); ++r) {
		for (auto n = info::mesh->elementsRegions[r]->elements->datatarray().begin(); n != info::mesh->elementsRegions[r]->elements->datatarray().end(); ++n) {
			*n = backpermutation[*n];
		}
		std::sort(info::mesh->elementsRegions[r]->elements->datatarray().begin(), info::mesh->elementsRegions[r]->elements->datatarray().end());
	}

	profiler::synccheckpoint("remap");
	profiler::syncend("permute_elements");
	eslog::checkpointln("MESH: ELEMENTS PERMUTED");
}

}
}
