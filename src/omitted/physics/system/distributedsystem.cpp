
#include "distributedsystem.h"
#include "builder/builder.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "basis/utilities/debugprint.h"
#include "basis/utilities/sysutils.h"
#include "wrappers/mpi/communication.h"

#include "math/matrix.csr.distributed.h"

#include <algorithm>
#include <fstream>

using namespace espreso;

void DistributedAssemblerData::print(const Builder *builder, const char* prefix, const char* suffix)
{
    if (builder->matrices & Builder::Request::K) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("K") + std::string(suffix)));
        os << K;
    }
    if (builder->matrices & Builder::Request::C) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("C") + std::string(suffix)));
        os << C;
    }
    if (builder->matrices & Builder::Request::C) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("CM") + std::string(suffix)));
        os << CM;
    }
    if (builder->matrices & Builder::Request::M) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("M") + std::string(suffix)));
        os << M;
    }
    if (builder->matrices & Builder::Request::R) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("R") + std::string(suffix)));
        os << R;
    }
    if (builder->matrices & Builder::Request::f) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("f") + std::string(suffix)));
        os << f;
    }
    if (builder->matrices & Builder::Request::BC) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("BC") + std::string(suffix)));
        os << BC;
    }
}

struct __Dirichlet__ {
    double value;
    esint row, column;
};

void DistributedSolverData::setDirichlet(const Builder *builder)
{
    if (builder->matrices & Builder::Request::KCM) {
        _KDirichletValues.clear();
    }
    if (!(builder->matrices & (Builder::Request::KCM | Builder::Request::f))) {
        return;
    }

    std::vector<__Dirichlet__> tosend;

    auto getrow = [&] (esint dof) {
        esint row = -1;
        if (K.distribution[info::mpi::rank] <= dof && dof < K.distribution[info::mpi::rank + 1]) {
            row = dof - K.distribution[info::mpi::rank] + K.nhalo;
        } else {
            auto inhalo = std::lower_bound(K.halo, K.halo + K.nhalo, dof);
            if (inhalo != K.halo + K.nhalo && *inhalo == dof) {
                row = inhalo - K.halo;
            }
        }
        return row;
    };

    size_t v = 0;
    for (esint i = 0; i < BC[0].nnz; ++i) {
        f[0].vals[BC[0].indices[i]] = BC[0].vals[i];
        esint col = 0;
        if (BC[0].indices[i] < K.nhalo) {
            col = K.halo[BC[0].indices[i]] + 1;
        } else {
            col = K.distribution[info::mpi::rank] + BC[0].indices[i] - K.nhalo + 1;
        }
        for (esint j = K.rows[BC[0].indices[i]]; j < K.rows[BC[0].indices[i] + 1]; j++) {
            if (K.cols[j - 1] == col) {
                K.vals[j - 1] = 1;
            } else {
                K.vals[j - 1] = 0;
                esint row = getrow(K.cols[j - 1] - 1);
                if (row != -1 && !std::binary_search(BC[0].indices, BC[0].indices + BC[0].nnz, row)) {
                    for (esint c = K.rows[row]; c < K.rows[row + 1]; c++) {
                        if (K.cols[c - 1] == col) {
                            double val;
                            if (v == _KDirichletValues.size()) {
                                val = K.vals[c - 1];
                                _KDirichletValues.push_back(val);
                            } else {
                                val = _KDirichletValues[v];
                            }
                            ++v;
                            if (row < K.nhalo) {
                                tosend.push_back({val * f[0].vals[BC[0].indices[i]], K.cols[j - 1] - 1, K.cols[c - 1] - 1});
                            }
                            f[0].vals[row] -= val * f[0].vals[BC[0].indices[i]];
                            K.vals[c - 1] = 0;
                        }
                    }
                } else {
                    // the column will be updated by a higher process
                }
            }
        }
    }

    std::sort(tosend.begin(), tosend.end(), [] (const __Dirichlet__ &i, const __Dirichlet__ &j) {
        if (i.row == j.row) {
            return i.column < j.column;
        }
        return i.row < j.row;
    });

    std::vector<std::vector<__Dirichlet__> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());
    for (size_t n = 0, i = 0; n < info::mesh->neighbors.size(); n++) {
        while (i < tosend.size() && tosend[i].row < K.distribution[info::mesh->neighbors[n] + 1]) {
            sBuffer[n].push_back(tosend[i++]);
        }
    }

    if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
        eslog::internalFailure("synchronize dirichlet data.\n");
    }

    for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
        for (size_t i = 0; i < rBuffer[n].size(); i++) {
            esint row = getrow(rBuffer[n][i].column);
            if (row == -1) {
                row = getrow(rBuffer[n][i].row);
                f[0].vals[row] -= rBuffer[n][i].value;
                esint c = std::lower_bound(K.cols + K.rows[row] - 1, K.cols + K.rows[row + 1] - 1, rBuffer[n][i].column + 1) - K.cols;
                K.vals[c] = 0;
            }
        }
    }
}

void DistributedSolverData::printData(const Builder *builder, const char* prefix)
{
    if (builder->matrices & Builder::Request::K) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("K")));
        os << K;
    }

    if (builder->matrices & Builder::Request::R) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("R")));
        os << R;
    }

    if (builder->matrices & Builder::Request::RBCf) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("f")));
        os << f[0];
    }

    if (builder->matrices & Builder::Request::BC) {
        std::ofstream os(utils::prepareFile(std::string(prefix), std::string("BC")));
        os << BC;
    }
}

void DistributedSolverData::printSolution(const Builder *builder, const char* prefix)
{
    std::ofstream os(utils::prepareFile(std::string(prefix), std::string("x")));
    os << x;
}




