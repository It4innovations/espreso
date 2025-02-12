
#include "distributed.composer.opt.h"

#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "physics/system/distributedsystem.h"
#include "physics/assembler/modules/module.opt.h"
#include "physics/kernels/solverdataprovider/provider.h"

using namespace espreso;

void DistributedComposerOpt::assemble(const Builder &builder)
{
    if (kernel) {
        DistributedComposer::assemble(builder);
        return;
    }
    if (builder.matrices == Builder::Request::NONE) {
        return;
    }

    eslog::startln("DISTRIBUTED ASSEMBLER: STARTED", "ASSEMBLER");

    clearMatrices(builder.matrices, _data);

    double start = eslog::time();

    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; t++) {
        for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
            bool omitLower = _data->K.type != MatrixType::REAL_UNSYMMETRIC;
            esint *Kperm = _KPermutation.data();
            esint *RHSperm = _RHSPermutation.data();

            for (esint ii = info::mesh->elements->eintervalsDistribution[d]; ii < info::mesh->elements->eintervalsDistribution[d + 1]; ++ii) {
                if (builder.matrices & Builder::Request::K) {
                    opt->updateStiffness(_data->K.vals, Kperm, ii);
                }
                if (builder.matrices & Builder::Request::M) {
//                    opt->updateStiffness(_data->M[d].vals, perm, ii);
                }
                if (builder.matrices & Builder::Request::f) {
//                    opt->updateStiffness(_data->M[d].vals, perm, ii);
                }
                Kperm += (info::mesh->elements->eintervals[ii].end - info::mesh->elements->eintervals[ii].begin) * getMatrixSize(esize(ii), omitLower);
                RHSperm += (info::mesh->elements->eintervals[ii].end - info::mesh->elements->eintervals[ii].begin) * esize(ii);
            }
            for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
                if (info::mesh->boundaryRegions[r]->dimension) {
                    for (esint ii = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; ii < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++ii) {
                        if (builder.matrices & Builder::Request::K) {
    //                        opt->updateStiffness(_data->K[d].vals, Kperm, ii);
                        }
                        if (builder.matrices & Builder::Request::M) {
        //                    opt->updateStiffness(_data->M[d].vals, perm, ii);
                        }
                        if (builder.matrices & Builder::Request::f) {
//                            for (esint n = 0; n < _data->f.nvectors; n++) {
//                                opt->updateRHS(_data->f[0].vals, RHSperm, r, ii);
//                            }
                        }
                        RHSperm += (info::mesh->boundaryRegions[r]->eintervals[ii].end - info::mesh->boundaryRegions[r]->eintervals[ii].begin) * bsize(r, ii);
                    }
                }
            }
        }
    }

    printf("INSERT: %f\n", eslog::time() - start);

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
        provider()->general->dirichletValues(values);
        for (size_t i = 0; i < values.size(); i++) {
            values[i] *= builder.internalForceReduction;
        }
        fillPermutedSparseData(_data->BC[0].vals, _dirichletMap, _dirichletPermutation, values);
        eslog::checkpointln("ASSEMBLER: DIRICHLET COMPUTED");
    }

    eslog::endln("ASSEMBLER: FINISHED");
}
