
#include "feti.composer.opt.h"

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
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "physics/system/fetisystem.h"
#include "physics/assembler/modules/module.opt.h"
#include "physics/kernels/solverdataprovider/provider.h"
#include "physics/system/builder/builder.h"
#include "math/vector.sparse.h"
#include "wrappers/bem/w.bem.h"
#include <cstddef>
#include <cmath>

using namespace espreso;

void FETIComposerOpt::assemble(const Builder &builder)
{
	if (kernel) {
		FETIComposer::assemble(builder);
		return;
	}
	if (builder.matrices == Builder::Request::NONE) {
		return;
	}

	eslog::startln("FETI ASSEMBLER: STARTED", "ASSEMBLER");

	clearMatrices(builder.matrices, _data);

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			bool omitLower = _data->K[d].type != MatrixType::REAL_UNSYMMETRIC;
			esint *Kperm = _KPermutation[d].data();
			esint *RHSperm = _RHSPermutation[d].data();

			for (esint ii = info::mesh->elements->eintervalsDistribution[d]; ii < info::mesh->elements->eintervalsDistribution[d + 1]; ++ii) {
				if (builder.matrices & Builder::Request::K) {
					opt->updateStiffness(_data->K[d].vals, Kperm, ii);
				}
				if (builder.matrices & Builder::Request::M) {
//					opt->updateStiffness(_data->M[d].vals, perm, ii);
				}
				if (builder.matrices & Builder::Request::f) {
//					opt->updateStiffness(_data->M[d].vals, perm, ii);
				}
				Kperm += (info::mesh->elements->eintervals[ii].end - info::mesh->elements->eintervals[ii].begin) * getMatrixSize(esize(ii), omitLower);
				RHSperm += (info::mesh->elements->eintervals[ii].end - info::mesh->elements->eintervals[ii].begin) * esize(ii);
			}
			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for (esint ii = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; ii < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++ii) {
						if (builder.matrices & Builder::Request::K) {
	//						opt->updateStiffness(_data->K[d].vals, Kperm, ii);
						}
						if (builder.matrices & Builder::Request::M) {
		//					opt->updateStiffness(_data->M[d].vals, perm, ii);
						}
						if (builder.matrices & Builder::Request::f) {
							opt->updateRHS(_data->f[0][d].vals, RHSperm, r, ii);
						}
						RHSperm += (info::mesh->boundaryRegions[r]->eintervals[ii].end - info::mesh->boundaryRegions[r]->eintervals[ii].begin) * bsize(r, ii);
					}
				}
			}
		}
	}

	if (builder.matrices & Builder::Request::BC) {
		std::vector<double> values(_dirichletMap.size());
		provider()->general->dirichletValues(values);

		for (size_t i = 0; i < values.size(); i++) {
			values[i] *= builder.internalForceReduction;
		}

		fillPermutedSparseData(_data->BC[0].vals, _dirichletMap, _dirichletPermutation, values);
		eslog::checkpointln("ASSEMBLER: DIRICHLET COMPUTED");
	}
	eslog::endln("ASSEMBLER: FINISHED");
}

