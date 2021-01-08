
#include "copy.h"

#include "basis/evaluator/evaluator.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "physics/assembler/parameter.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

void AverageElementsNodesToNodes::now()
{
	int version = *std::max_element(from.version.begin(), from.version.begin());
	Communication::allReduce(&version, NULL, 1, MPI_INT, MPI_MAX);
	if (to.version < version || !to.version) {
		if (Operator::print) printf("\tOP::AverageElementsNodesToNodes\n");
		to.version = version;
		std::fill(to.data.begin(), to.data.end(), 0);
		auto procNodes = info::mesh->elements->procNodes->begin();
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			InputParameterIterator input(from, i, to.dimension);
			for (esint e = info::mesh->elements->eintervals[i].begin; e < info::mesh->elements->eintervals[i].end; ++e, ++procNodes) {
				for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++input) {
					for (int d = 0; d < to.dimension; ++d) {
						to.data[*n * to.dimension + d] += input.data[d];
					}
				}
			}
		}

		auto &nelements = info::mesh->nodes->elements->boundarytarray();
		for (size_t i = 0; i < to.data.size() / to.dimension; i++) {
			for (int d = 0; d < to.dimension; d++) {
				to.data[i * to.dimension + d] /= nelements[i + 1] - nelements[i];
			}
		}

		std::vector<std::vector<double> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());

		auto nranks = info::mesh->nodes->ranks->begin();
		for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
			if (nranks->size() > 1) {
				esint noffset = 0;
				for (auto r = nranks->begin(); r != nranks->end(); ++r) {
					if (*r != info::mpi::rank) {
						while (info::mesh->neighbors[noffset] < *r) {
							++noffset;
						}
						for (int d = 0; d < to.dimension; d++) {
							sBuffer[noffset].push_back(to.data[n * to.dimension + d]);
						}
					}
				}
			}
		}

		for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
			rBuffer[n].resize(sBuffer[n].size());
		}

		if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
			eslog::internalFailure("exchange diagonal values.\n");
		}

		nranks = info::mesh->nodes->ranks->begin();
		std::vector<esint> nindex(info::mesh->neighbors.size());
		for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
			if (nranks->size() > 1) {
				esint noffset = 0;
				for (auto r = nranks->begin(); r != nranks->end(); ++r) {
					if (*r != info::mpi::rank) {
						while (info::mesh->neighbors[noffset] < *r) {
							++noffset;
						}
						for (int d = 0; d < to.dimension; d++) {
							to.data[n * to.dimension + d] += rBuffer[noffset][nindex[noffset]++];
						}
					}
				}
			}
		}
	} else {
		printf("\tOP::SKIPPED::AverageElementsNodesToNodes\n");
	}
}

void CopyNodesToElementsNodes::now()
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; t++) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			for (esint ii = info::mesh->elements->eintervalsDistribution[d]; ii < info::mesh->elements->eintervalsDistribution[d + 1]; ++ii) {
				if (to.version[ii] < from.version || !to.version[ii]) {
					to.version[ii] = from.version;
					if (Operator::print) printf("\tOP::CopyNodesToElementsNodes::%d\n", ii);
					auto i = (to.data->begin() + ii)->data();
					auto procNodes = info::mesh->elements->procNodes->begin() + info::mesh->elements->eintervals[ii].begin;
					for (esint e = info::mesh->elements->eintervals[ii].begin; e < info::mesh->elements->eintervals[ii].end; ++e, ++procNodes) {
						for (auto n = procNodes->begin(); n != procNodes->end(); ++n) {
							for (int dim = 0; dim < from.dimension; ++dim, ++i) {
								*i = from.data[*n * from.dimension + dim];
							}
						}
					}
				} else {
					if (Operator::print) printf("\tOP::CopyNodesToElementsNodes::%d::SKIPPED\n", ii);
				}
			}
		}
	}
}

void CopyNodesToBoundaryNodes::now()
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; t++) {
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; ++d) {
					for (esint ii = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; ii < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++ii) {
						if (to.regions[r].version[ii] < from.version || !to.regions[r].version[ii]) {
							to.regions[r].version[ii] = from.version;
							if (Operator::print) printf("\tOP::CopyNodesToBoundaryNodes::%lu::%d\n", r, ii);

							auto i = (to.regions[r].data->begin() + ii)->data();
							auto procNodes = info::mesh->boundaryRegions[r]->procNodes->begin() + info::mesh->boundaryRegions[r]->eintervals[ii].begin;
							for (auto e = info::mesh->boundaryRegions[r]->eintervals[ii].begin; e != info::mesh->boundaryRegions[r]->eintervals[ii].end; ++e, ++procNodes) {
								for (auto n = procNodes->begin(); n != procNodes->end(); ++n) {
									for (int dim = 0; dim < from.dimension; ++dim, ++i) {
										*i = from.data[*n * from.dimension + dim];
									}
								}
							}
						}  else {
							if (Operator::print) printf("\tOP::CopyNodesToBoundaryNodes::%lu::%d::SKIPPED\n", r, ii);
						}
					}
				}
			}
		}
	}
}

void CopyBoundaryRegionsSettingToNodes::now()
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		for (auto it = from.begin(); it != from.end(); ++it) {
			serializededata<esint, esint>* nodes = info::mesh->bregion(it->first)->nodes;
			it->second.evaluator->evalSelectedDense(nodes->datatarray().size(t), nodes->datatarray().begin(t), it->second.evaluator->params, to.data.data());
		}
	}
	++to.version;
}

