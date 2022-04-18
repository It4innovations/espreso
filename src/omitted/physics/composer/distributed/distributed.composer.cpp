
#include "distributed.composer.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/domainstore.h"
#include "physics/system/distributedsystem.h"
#include "physics/kernels/kernel.h"
#include "physics/kernels/solverdataprovider/provider.h"
#include "physics/system/builder/builder.h"
#include "math/vector.sparse.h"

#include <cstddef>
#include <cmath>

using namespace espreso;

DistributedComposer::DistributedComposer(Kernel *kernel, ModuleOpt *opt, DistributedAssemblerData *data)
: Composer(kernel, opt), _data(data),
  _DOFMap(NULL)
{

}

DistributedComposer::~DistributedComposer()
{
	if (_DOFMap != NULL) { delete _DOFMap; }
}

void DistributedComposer::assemble(const Builder &builder)
{
	if (builder.matrices == Builder::Request::NONE) {
		return;
	}

	eslog::startln("DISTRIBUTED ASSEMBLER: STARTED", "ASSEMBLER");

	int invalid = 0;

	clearMatrices(builder.matrices, _data);

	double insertTime = 0, assembleTime = 0;

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; t++) {
		size_t KIndex = _tKOffsets[t], RHSIndex = _tRHSOffsets[t];
		double KReduction = builder.timeIntegrationConstantK, RHSReduction = builder.internalForceReduction;
		Kernel::InstanceFiller filler(solutions());

		double prev, tinsertTime = 0, tassembleTime = 0;
		filler.insert = [&] () {
			tassembleTime += eslog::time() - prev;
			prev = eslog::time();
			for (int r = 0; r < filler.DOFs; ++r, ++RHSIndex) {
				if (filler.insertF) {
					for (esint n = 0; n < _data->f.nvectors; n++) {
						#pragma omp atomic
						_data->f[n].vals[_RHSPermutation[RHSIndex]] += RHSReduction * filler.Fe[n].vals[r];
					}
				}
				if (filler.insertR) {
					for (esint n = 0; n < _data->R.nvectors; n++) {
						#pragma omp atomic
						_data->R[n].vals[_RHSPermutation[RHSIndex]] += filler.Re[n].vals[r];
					}
				}

				for (int c = 0; c < filler.DOFs; ++c, ++KIndex) {
					if (filler.insertK) {
						#pragma omp atomic
						_data->K.vals[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
					}
					if (filler.insertM) {
						#pragma omp atomic
						_data->M.vals[_KPermutation[KIndex]] += filler.Me(r, c);
					}
					if (filler.insertC) {
						#pragma omp atomic
						_data->C.vals[_KPermutation[KIndex]] += filler.Ce(r, c);
						#pragma omp atomic
						_data->CM.vals[_KPermutation[KIndex]] += filler.CMe(r, c);
					}
				}
			}
			tinsertTime += eslog::time() - prev;
			prev = eslog::time();
		};

		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint ii = info::mesh->elements->eintervalsDistribution[d]; ii < info::mesh->elements->eintervalsDistribution[d + 1]; ++ii) {
				filler.interval = ii;
				prev = eslog::time();
				kernel->processElements(builder, filler);
			}
		}

		if (t == 0) {
			eslog::checkpointln("ASSEMBLER: ELEMENTS PROCESSED");
		}

		KReduction = builder.internalForceReduction;
		filler.insertM = filler.insertC = filler.insertR = false;

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->dimension) {
                                filler.begin = info::mesh->boundaryRegions[r]->epointers->datatarray().distribution()[t];
                                filler.end = info::mesh->boundaryRegions[r]->epointers->datatarray().distribution()[t + 1];
				prev = eslog::time();
				kernel->processBoundary(builder, r, filler);
			}
		}

		if (t == 0) {
			eslog::checkpointln("ASSEMBLER: BOUNDARY REGIONS PROCESSED");
		}

		esint noffset = info::mesh->nodes->size - info::mesh->nodes->uniqInfo.size;
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (!info::mesh->boundaryRegions[r]->dimension) {
				// process only nodes held by this process
				esint roffset = std::lower_bound(
						info::mesh->boundaryRegions[r]->nodes->datatarray().begin(),
						info::mesh->boundaryRegions[r]->nodes->datatarray().end(),
						noffset) - info::mesh->boundaryRegions[r]->nodes->datatarray().begin();
				esint begin = info::mesh->boundaryRegions[r]->nodes->datatarray().distribution()[t];
				esint end = info::mesh->boundaryRegions[r]->nodes->datatarray().distribution()[t + 1];
				if (roffset < end) {
					filler.begin = std::max(roffset, begin);
					filler.end = end;
					prev = eslog::time();
					kernel->processBoundary(builder, r, filler);
				}
			}
		}

		if (t == 0) {
			eslog::checkpointln("ASSEMBLER: NODES REGIONS PROCESSED");
		}

		#pragma omp atomic
		invalid += filler.invalid;
		#pragma omp atomic
		assembleTime += tassembleTime;
		#pragma omp atomic
		insertTime += tinsertTime;
	}

	eslog::checkpoint("ASSEMBLER: VALUES FILLED");
	eslog::param("ASSEMBLE", assembleTime);
	eslog::param("INSERT", insertTime);
	eslog::ln();

	int ginvalid;
	Communication::reduce(&invalid, &ginvalid, 1, MPI_INT, MPI_SUM, 0);
	if (info::mpi::rank == 0 && ginvalid) {
		eslog::warning("ESPRESO mesh warning: %d elements with negative determinant.\n", ginvalid);
	}

	if (builder.matrices & Builder::Request::K) {
		_data->K.gatherFromUpper();
	}
	if (builder.matrices & Builder::Request::M) {
		_data->M.gatherFromUpper();
	}
	if (builder.matrices & Builder::Request::C) {
		_data->C.gatherFromUpper();
	}
	if (builder.matrices & Builder::Request::R) {
		_data->R.gatherFromUpper();
	}
	if (builder.matrices & Builder::Request::f) {
		_data->f.gatherFromUpper();
	}

	eslog::checkpointln("ASSEMBLER: DATA GATHERED");

	if (builder.matrices & Builder::Request::BC) {
		std::vector<double> values(_dirichletMap.size());
		kernel->solverDataProvider->general->dirichletValues(values);
		for (size_t i = 0; i < values.size(); i++) {
			values[i] *= builder.internalForceReduction;
		}
		fillPermutedSparseData(_data->BC[0].vals, _dirichletMap, _dirichletPermutation, values);
		eslog::checkpointln("ASSEMBLER: DIRICHLET COMPUTED");
	}

	eslog::endln("ASSEMBLER: FINISHED");
}


