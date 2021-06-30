
#include "feti.composer.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "physics/system/fetisystem.h"
#include "physics/kernels/kernel.h"
#include "physics/kernels/solverdataprovider/provider.h"
#include "physics/system/builder/builder.h"
#include "math/vector.sparse.h"
#include "wrappers/bem/w.bem.h"
#include <cstddef>
#include <cmath>

#include "basis/utilities/print.h"

using namespace espreso;

FETIComposer::FETIComposer(const FETIConfiguration &configuration, Kernel *kernel, ModuleOpt *opt, FETIAssemblerData *data)
: Composer(kernel, opt), _configuration(configuration), _data(data), _DOFMap(NULL)
{
//	auto isBEMDomain = [&](esint domain) {
//		auto eregions = (info::mesh->elements->regions->begin() + info::mesh->domains->elements[domain])->begin();
//		for (int byte = 0; byte < info::mesh->elements->regionMaskSize; ++byte) {
//			for (size_t bit = 0; bit < sizeof(esint); bit++) {
//				if (eregions[byte] & 1 << bit) {
//					auto region = info::mesh->elementsRegions[byte * sizeof(esint) + bit];
//					if (assembler->configuration().discretization.find(region->name) != assembler->configuration().discretization.end()) {
//						return true;
//					}
//				}
//			}
//		}
//		return false;
//	};

	_BEMDomain.resize(info::mesh->domains->size);
	if (BEM4I::isLinked()) {
//		for (esint d = 0; d < info::mesh->domains->size; ++d) {
//			_BEMDomain[d] = isBEMDomain(d);
//		}
	}
}

FETIComposer::~FETIComposer()
{
	if (_DOFMap) { delete _DOFMap; }
}

void FETIComposer::assemble(const Builder &builder)
{
	if (builder.matrices == Builder::Request::NONE) {
		return;
	}

	eslog::startln("FETI ASSEMBLER: STARTED", "ASSEMBLER");

	clearMatrices(builder.matrices, _data);
	int invalid = 0;
	double insertTime = 0, assembleTime = 0;

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			size_t KIndex = 0, RHSIndex = 0;
			double KReduction = builder.timeIntegrationConstantK, RHSReduction = builder.internalForceReduction;
			Kernel::InstanceFiller filler(solutions());
			filler.K = _data->K[d].vals;
			filler.M = _data->M[d].vals;
			filler.C = _data->C[d].vals;
			filler.R = _data->R[0][d].vals;
			filler.F = _data->f[0][d].vals;
			filler.offset = _KPermutation[d].data();

			double prev, tinsertTime = 0, tassembleTime = 0;
			switch (kernel->solverDataProvider->feti->getMatrixType(d)) {
			case MatrixType::REAL_UNSYMMETRIC:
				filler.insert = [&] () {
					tassembleTime += eslog::time() - prev;
					prev = eslog::time();
					for (int r = 0; r < filler.DOFs; ++r, ++RHSIndex) {
						if (filler.insertF) {
							for (esint n = 0; n < _data->f.nvectors; n++) {
								_data->f[n][d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.Fe[n][r];
							}
						}
						if (filler.insertR) {
							for (esint n = 0; n < _data->R.nvectors; n++) {
								_data->R[n][d][_RHSPermutation[d][RHSIndex]] += filler.Re[n][r];
							}
						}
						for (int c = 0; c < filler.DOFs; ++c, ++KIndex) {
							if (filler.insertK) {
								_data->K[d].vals[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
							}
							if (filler.insertM) {
								_data->M[d].vals[_KPermutation[d][KIndex]] += filler.Me(r % filler.Me.nrows, c % filler.Me.nrows);
							}
							if (filler.insertC) {
								_data->C[d].vals[_KPermutation[d][KIndex]] += filler.Ce(r, c);
								_data->CM[d].vals[_KPermutation[d][KIndex]] += filler.CMe(r, c);
							}
						}
					}
					tinsertTime += eslog::time() - prev;
					prev = eslog::time();
				}; break;
			default:
				filler.insert = [&] () {
					tassembleTime += eslog::time() - prev;
					prev = eslog::time();
					for (int r = 0; r < filler.DOFs; ++r, ++RHSIndex) {
						if (filler.insertF) {
							for (esint n = 0; n < _data->f.nvectors; n++) {
								_data->f[n][d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.Fe[n][r];
							}
						}
						if (filler.insertR) {
							for (esint n = 0; n < _data->R.nvectors; n++) {
								_data->R[n][d][_RHSPermutation[d][RHSIndex]] += filler.Re[n][r];
							}
						}

						for (int c = r; c < filler.DOFs; ++c, ++KIndex) {
							if (filler.insertK) {
								_data->K[d].vals[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
							}
							if (filler.insertM) {
								_data->M[d].vals[_KPermutation[d][KIndex]] += filler.Me(r % filler.Me.nrows, c % filler.Me.nrows);
							}
							if (filler.insertC) {
								_data->C[d].vals[_KPermutation[d][KIndex]] += filler.Ce(r, c);
								_data->CM[d].vals[_KPermutation[d][KIndex]] += filler.CMe(r, c);
							}
						}
					}
					tinsertTime += eslog::time() - prev;
					prev = eslog::time();
				}; break;
			}

			for (esint ii = info::mesh->elements->eintervalsDistribution[d]; ii < info::mesh->elements->eintervalsDistribution[d + 1]; ++ii) {
				filler.interval = ii;
				if (_BEMDomain[d]) {
		//			assembler->processBEMdomain(d, _data->K[d].vals);
				} else {
					prev = eslog::time();
					kernel->processElements(builder, filler);
				}
			}

			KReduction = builder.internalForceReduction;
			filler.insertM = filler.insertC = filler.insertR = false;

			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					if (info::mesh->boundaryRegions[r]->eintervalsDistribution[d] < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]) {
						filler.begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[d]].begin;
						filler.end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
						prev = eslog::time();
						kernel->processBoundary(builder, r, filler);
					}
				}
			}
			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
				if (!info::mesh->boundaryRegions[r]->dimension) {
					esint prev = 0, i = 0;
					auto dmap = info::mesh->nodes->domains->begin();
					for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); prev = *n++, ++i) {
						dmap += *n - prev;
						if (dmap->at(0) == (esint)d + info::mesh->domains->offset) {
							filler.begin = i;
							filler.end = i + 1;
							prev = eslog::time();
							kernel->processBoundary(builder, r, filler);
						}
					}
				}
			}

			#pragma omp atomic
			invalid += filler.invalid;
			#pragma omp atomic
			assembleTime += tassembleTime;
			#pragma omp atomic
			insertTime += tinsertTime;
		}
	}

	printf("ASSEMBLE: %.3f, INSERT: %.3f\n", assembleTime, insertTime);

	eslog::checkpoint("ASSEMBLER: VALUES FILLED");
	eslog::param("ASSEMBLE", assembleTime);
	eslog::param("INSERT", insertTime);
	eslog::ln();

	int ginvalid;
	Communication::reduce(&invalid, &ginvalid, 1, MPI_INT, MPI_SUM, 0);
	if (info::mpi::rank == 0 && ginvalid) {
		eslog::warning("ESPRESO mesh warning: %d elements with negative determinant.\n", ginvalid);
	}

//	synchronize(builder);

	if (builder.matrices & Builder::Request::BC) {
		std::vector<double> values(_dirichletMap.size());
		kernel->solverDataProvider->general->dirichletValues(values);

		for (size_t i = 0; i < values.size(); i++) {
			values[i] *= builder.internalForceReduction;
		}

		fillPermutedSparseData(_data->BC[0].vals, _dirichletMap, _dirichletPermutation, values);
		eslog::checkpointln("ASSEMBLER: DIRICHLET COMPUTED");
	}
//
//	if (builder.matrices & Builder::Request::BC) {
//		std::vector<double> values(_dirichletMap.size());
//		assembler->dirichletValues(values);
//
//		for (size_t i = 0; i < values.size(); i++) {
//			values[i] *= builder.internalForceReduction;
//		}
//
//		auto map = _data->BC.holder()->sparsemap->begin();
//		for (size_t i = 0, j = 0; i < _dirichletMap.size(); i = j, ++map) {
//			for (auto di = map->begin(); di != map->end(); ++di) {
//				if (_data->BC.holder()->ismy(di->domain)) {
//					_data->BC.holder()->at(di->domain - _data->BC.holder()->doffset)->vals[di->index] = 0;
//				}
//			}
//			while (j < _dirichletMap.size() && _dirichletMap[j] == _dirichletMap[i]) {
//				for (auto di = map->begin(); di != map->end(); ++di) {
//					if (_data->BC.holder()->ismy(di->domain)) {
//						_data->BC.holder()->at(di->domain - _data->BC.holder()->doffset)->vals[di->index] += values[_dirichletPermutation[j]];
//					}
//				}
//				++j;
//			}
//			for (auto di = map->begin(); di != map->end(); ++di) {
//				if (_data->BC.holder()->ismy(di->domain)) {
//					_data->BC.holder()->at(di->domain - _data->BC.holder()->doffset)->vals[di->index] /= j - i;
//				}
//			}
//		}
//	}
	eslog::endln("ASSEMBLER: FINISHED");
}


