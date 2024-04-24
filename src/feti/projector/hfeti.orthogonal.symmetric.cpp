
#include "hfeti.orthogonal.symmetric.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"

#include <vector>
#include <unordered_map>

namespace espreso {

template<typename T>
HFETIOrthogonalSymmetric<T>::HFETIOrthogonalSymmetric(FETI<T> &feti)
: Projector<T>(feti)
{
    iGGtGx.resize(feti.sinfo.R1size);

    cinfo.reserve(1);
    cinfo.push_back(ClusterInfo(info::mpi::rank, feti.sinfo.R1offset, feti.sinfo.R1size));

    _computeDualGraph();
    _setG();
    _setGGt();
}

template<typename T>
HFETIOrthogonalSymmetric<T>::~HFETIOrthogonalSymmetric()
{

}

template<typename T>
void HFETIOrthogonalSymmetric<T>::info()
{
    int nnz = 2 * (GGt.nnz - GGt.nrows) + GGt.nrows;

    eslog::info(" = ORTHOGONAL PROJECTOR PROPERTIES                                                           = \n");
    eslog::info(" =   GGT ROWS                                                                      %9d = \n", GGt.nrows);
    eslog::info(" =   GGT NNZ                                                                       %9d = \n", nnz);
//    eslog::info(" =   GGT FACTORS NNZ                                                               %9d = \n", GGtSolver.nnzL);
    if (feti.configuration.exhaustive_info) {
        // PPt = eye
    }
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::update(const step::Step &step)
{
    Vector_Dense<T> _e;
    _e.size = feti.sinfo.R1size;
    _e.vals = e.vals + feti.sinfo.R1offset;
    math::set(e, T{0});
    for (size_t d = 0; d < feti.R1.size(); ++d) {
        math::blas::apply(_e, T{-1}, feti.R1[d], T{1}, feti.f[d]);
    }
    e.synchronize();
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [e]");

    if (feti.updated.K || feti.updated.B) {
        _updateG();
        _updateGGt();
    }
    _print(step);
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    if (GGt.nrows) {
        x.copyToWithoutHalo(y);
        _applyG(x, Gx);
        _applyInvGGt(Gx, iGGtGx);
        _applyGt(iGGtGx, T{-1}, y);
    } else {
        math::copy(y, x);
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply_e(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
    if (GGt.nrows) {
        math::set(y, T{0});
        _applyInvGGt(x, iGGtGx);
        _applyGt(iGGtGx, T{1}, y);
    } else {
        math::set(y, T{0});
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply_R(const Vector_Kernel<T> &x,  std::vector<Vector_Dense<T> > &y)
{
    if (GGt.nrows) {
        _applyR(x, y);
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < y.size(); ++d) {
            math::set(y[d], T{0});
        }
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply_Ra(const Vector_Dual<T> &x,  std::vector<Vector_Dense<T> > &y)
{
    if (GGt.nrows) {
        _applyG(x, Gx);
        _applyInvGGt(Gx, iGGtGx);
        _applyR(iGGtGx, y);
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < y.size(); ++d) {
            math::set(y[d], T{0});
        }
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply_invU(const Vector_Kernel<T> &x, Vector_Kernel<T> &y)
{
    if (GGt.nrows) {
        _applyInvU(x, y);
    } else {
        math::set(y, T{0});
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply_invL(const Vector_Kernel<T> &x, Vector_Kernel<T> &y)
{
    if (GGt.nrows) {
        _applyInvL(x, y);
    } else {
        math::set(y, T{0});
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply_GtinvU(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
    if (GGt.nrows) {
        math::set(y, T{0});
        _applyInvU(x, iGGtGx);
        _applyGt(iGGtGx, T{1}, y);
    } else {
        math::set(y, T{0});
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::apply_invLG(const Vector_Dual<T> &x, Vector_Kernel<T> &y)
{
    if (GGt.nrows) {
        _applyG(x, Gx);
        _applyInvL(Gx, y);
    } else {
        math::copy(y, x);
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out)
{
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t r = Vector_Kernel<T>::distribution[t]; r < Vector_Kernel<T>::distribution[t + 1]; ++r) {
            out.vals[r + Vector_Kernel<T>::offset] = T{0};
            for (int c = 0; c < G.ncols; ++c) {
                out.vals[r + Vector_Kernel<T>::offset] += G.vals[r * G.ncols + c] * in.vals[c];
            }
        }
    }
    out.synchronize();
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_applyInvGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        Matrix_Dense<T> a;
        Vector_Dense<T> y;
        a.ncols = invGGt.ncols;
        a.nrows = y.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];

        a.vals = invGGt.vals + invGGt.ncols * Vector_Kernel<T>::distribution[t];
        y.vals = out.vals + Vector_Kernel<T>::distribution[t];

        math::blas::apply(y, T{1}, a, T{0}, in);
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_applyInvL(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
    eslog::error("cannot apply inv(L). Use different projector.\n");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_applyInvU(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
    eslog::error("cannot apply inv(U). Use different projector.\n");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out)
{
    for (int r = 0; r < G.nrows; ++r) {
        for (int c = 0; c < G.ncols; ++c) {
            out.vals[c] += alpha * G.vals[r * G.ncols + c] * in.vals[r];
        }
    }
    out.synchronize();
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_applyR(const Vector_Dense<T> &in, std::vector<Vector_Dense<T> > &out)
{
    #pragma omp parallel for
    for (size_t d = 0; d < out.size(); ++d) {
        Vector_Dense<T> y;
        y.size = cinfo[0].kernels;
        y.vals = in.vals + cinfo[0].koffset - feti.sinfo.R1offset;

        math::blas::applyT(out[d], T{1}, feti.R1[d], T{0}, y);
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_computeDualGraph()
{
    std::vector<std::vector<ClusterInfo> > sBuffer(feti.decomposition->neighbors.size()), rBuffer(feti.decomposition->neighbors.size());
    for (size_t n = 0; n < feti.decomposition->neighbors.size(); ++n) {
        sBuffer[n].push_back(cinfo.front());
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, feti.decomposition->neighbors)) {
        eslog::error("cannot exchange dual graph info\n");
    }

    dualGraph.resize(1);
    size_t n = 0;
    for (; n < feti.decomposition->neighbors.size() && feti.decomposition->neighbors[n] < info::mpi::rank; ++n) {
        dualGraph[0].push_back(ClusterInfo(feti.decomposition->neighbors[n], rBuffer[n][0].koffset, rBuffer[n][0].kernels));
    }
    dualGraph[0].push_back(ClusterInfo(info::mpi::rank, feti.sinfo.R1offset, feti.sinfo.R1size));
    for (; n < feti.decomposition->neighbors.size(); ++n) {
        dualGraph[0].push_back(ClusterInfo(feti.decomposition->neighbors[n], rBuffer[n][0].koffset, rBuffer[n][0].kernels));
    }

    neighInfo.resize(feti.decomposition->neighbors.size());
    for (size_t i = 0, offset = 0; i < feti.lambdas.cmap.size(); ) {
        int lsize = feti.lambdas.cmap[i];
        int domains = feti.lambdas.cmap[i + 1];
        for (int d = 0, last = -1; d < domains; ++d) {
            int n = feti.decomposition->noffset(feti.lambdas.cmap[i + 2 + d]);
            if (!feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d]) && last < n) {
                if (feti.lambdas.cmap[i + 2 + d] < feti.decomposition->dbegin) {
                    neighInfo[n] = dualGraph[0][n];
                    neighInfo[n].cindices.push_back({ (int)offset, lsize });
                    neighInfo[n].ncols += lsize;
                }
                if (feti.decomposition->dend <= feti.lambdas.cmap[i + 2 + d]) {
                    neighInfo[n] = dualGraph[0][n + 1];
                    neighInfo[n].cindices.push_back({ (int)offset, lsize });
                    neighInfo[n].ncols += lsize;
                }
                last = n;
            }
        }
        i += feti.lambdas.cmap[i + 1] + 2;
        offset += lsize;
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_setG()
{
    G.resize(feti.sinfo.R1size, feti.lambdas.size);

    int Gtrows = 0, Gtnnz = 0;
    for (size_t n = 0; n < feti.decomposition->neighbors.size(); ++n) {
        if (info::mpi::rank < feti.decomposition->neighbors[n]) {
            Gtrows += neighInfo[n].kernels;
            Gtnnz  += neighInfo[n].kernels * neighInfo[n].ncols;
        }
    }

    Gt.resize(Gtrows, feti.lambdas.size, Gtnnz);
    Gt.rows[0] = 0;
    for (size_t n = 0, ri = 0; n < feti.decomposition->neighbors.size(); ++n) {
        if (info::mpi::rank < feti.decomposition->neighbors[n]) {
            neighInfo[n].koffset = ri;
            for (int kr = 0; kr < neighInfo[n].kernels; ++kr, ++ri) {
                Gt.rows[ri + 1] = Gt.rows[ri] + neighInfo[n].ncols;
                for (size_t ci = 0, c = 0; ci < neighInfo[n].cindices.size(); ++ci) {
                    for (int cc = 0; cc < neighInfo[n].cindices[ci].count; ++cc, ++c) {
                        Gt.cols[Gt.rows[ri] + c] = neighInfo[n].cindices[ci].offset + cc;
                    }
                }
            }
        }
    }

    eslog::checkpointln("FETI: SET G");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_updateG()
{
    math::set(G, T{0});
    for (size_t d = 0; d < feti.K.size(); ++d) {
        for (int kr = 0; kr < feti.sinfo.R1size; ++kr) {
            for (int c = 0; c < feti.B1[d].nrows; ++c) {
                for (int i = feti.B1[d].rows[c]; i < feti.B1[d].rows[c + 1]; ++i) {
                    G.vals[kr * G.ncols + feti.D2C[d][c]] -= feti.R1[d].vals[feti.R1[d].ncols * kr + feti.B1[d].cols[i]] * feti.B1[d].vals[i];
                }
            }
        }
    }

    const FETIDecomposition *decomposition = feti.decomposition;
    std::vector<std::vector<T> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
    for (size_t n = 0; n < decomposition->neighbors.size(); ++n) {
        if (feti.decomposition->neighbors[n] < info::mpi::rank) {
            for (int kr = 0; kr < neighInfo[n].kernels; ++kr) {
                for (size_t cc = 0; cc < neighInfo[n].cindices.size(); ++cc) {
                    for (int c = 0; c < neighInfo[n].cindices[cc].count; ++c) {
                        sBuffer[n].push_back(G.vals[kr * G.ncols + neighInfo[n].cindices[cc].offset + c]);
                    }
                }
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
        eslog::error("cannot exchange neighbor's G\n");
    }

    for (size_t n = 0, offset = Gt.rows[0]; n < rBuffer.size(); ++n) {
        std::copy(rBuffer[n].begin(), rBuffer[n].end(), Gt.vals + offset);
        offset += rBuffer[n].size();
    }
    eslog::checkpointln("FETI: UPDATE G");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_setGGt()
{
    const int IDX = Indexing::CSR;
    GGtDataOffset = 0;
    for (size_t i = 0; i < dualGraph[0].size(); ++i) {
        for (int kr = 0; kr < cinfo[0].kernels; ++kr) {
            for (int kc = 0; kc < dualGraph[0][i].kernels; ++kc) {
                if (cinfo[0].koffset + kr <= dualGraph[0][i].koffset + kc) {
                    ++GGtDataOffset;
                }
            }
        }
    }
    GGtDataSize = GGtDataOffset;
    GGtNnz = Communication::exscan(GGtDataOffset);

    GGt.resize(feti.sinfo.R1totalSize, feti.sinfo.R1totalSize, GGtNnz);
    GGt.shape = Matrix_Shape::UPPER;
    GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    GGt.rows[0] = IDX;
    GGt.rows[feti.sinfo.R1offset] = GGtDataOffset + IDX;

    for (int kr = 0; kr < feti.sinfo.R1size; ++kr) {
        GGt.rows[feti.sinfo.R1offset + kr + 1] = GGt.rows[feti.sinfo.R1offset + kr];
        for (size_t i = 0, c = GGt.rows[feti.sinfo.R1offset + kr] - IDX; i < dualGraph[0].size(); ++i) {
            for (int kc = 0; kc < dualGraph[0][i].kernels; ++kc) {
                if (feti.sinfo.R1offset + kr <= dualGraph[0][i].koffset + kc) {
                    GGt.cols[c++] = dualGraph[0][i].koffset + kc + IDX;
                    ++GGt.rows[feti.sinfo.R1offset + kr + 1];
                }
            }
        }
    }

    if (!Communication::allGatherInplace(GGt.rows, feti.sinfo.R1offset + 1, G.nrows)) {
        eslog::error("cannot gather GGt rows.\n");
    }
    if (!Communication::allGatherInplace(GGt.cols, GGtDataOffset, GGtDataSize)) {
        eslog::error("cannot gather GGt cols.\n");
    }

    invGGt.resize(G.nrows, feti.sinfo.R1totalSize);
    eslog::checkpointln("FETI: SET GGT");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_updateGGt()
{
    const int IDX = Indexing::CSR;
    for (int kr = 0; kr < feti.sinfo.R1size; ++kr) {
        for (size_t i = 0, c = GGt.rows[cinfo[0].koffset + kr] - IDX; i < dualGraph[0].size(); ++i) {
            for (int kc = 0; kc < dualGraph[0][i].kernels; ++kc) {
                if (feti.sinfo.R1offset + kr <= dualGraph[0][i].koffset + kc) {
                    GGt.vals[c] = 0;
                    if (dualGraph[0][i].koffset - feti.sinfo.R1offset < G.nrows) { // local
                        Vector_Dense<T> r1; r1.size = G.ncols; r1.vals = G.vals + kr * G.ncols;
                        Vector_Dense<T> r2; r2.size = G.ncols; r2.vals = G.vals + kc * G.ncols;
                        math::blas::multiply(T{1}, r1, r2, T{0}, GGt.vals[c]);
                    } else {
                        int koffset = neighInfo[i - 1].koffset; // (i - i) since there is local cluster in dualGraph
                        for (int kk = Gt.rows[koffset + kc]; kk < Gt.rows[koffset + kc + 1]; ++kk) {
                            GGt.vals[c] += Gt.vals[kk] * G.vals[kr * G.ncols + Gt.cols[kk]];
                        }
                    }
                    ++c;
                }
            }
        }
    }

    if (!Communication::allGatherInplace(GGt.vals, GGtDataOffset, GGtDataSize)) {
        eslog::error("cannot gather GGt vals.\n");
    }
    eslog::checkpointln("FETI: GATHER GGT VALUES");

    if (GGt.nrows) {
        DirectSparseSolver<T> GGtSolver;
        GGtSolver.commit(GGt);
        GGtSolver.symbolicFactorization();
        GGtSolver.numericalFactorization();
        eslog::checkpointln("FETI: GGT FACTORIZATION");

        Matrix_Dense<T> eye;
        eye.resize(G.nrows, feti.sinfo.R1totalSize);
        math::set(eye, T{});
        for (int r = 0; r < G.nrows; ++r) {
            eye.vals[r * feti.sinfo.R1totalSize + feti.sinfo.R1offset + r] = T{1};
        }
        GGtSolver.solve(eye, invGGt);
    }
    eslog::checkpointln("FETI: COMPUTE GGT INVERSE");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_print(const step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/projector/{G, Gt, e, GGt, invGGt}\n");
        math::store(G, utils::filename(utils::debugDirectory(step) + "/feti/projector", "G").c_str());
        math::store(Gt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "Gt").c_str());
        math::store(e, utils::filename(utils::debugDirectory(step) + "/feti/projector", "e").c_str());
        math::store(GGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "GGt").c_str());
        math::store(invGGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invGGt").c_str());
    }
}

template struct HFETIOrthogonalSymmetric<double>;
template struct HFETIOrthogonalSymmetric<std::complex<double> >;

}
