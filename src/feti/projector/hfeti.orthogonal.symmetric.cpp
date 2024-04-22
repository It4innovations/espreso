
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

    domainOffset = feti.decomposition->dbegin;

    dinfo.reserve(feti.R1.size());
    for (size_t d = 0, koffset = feti.sinfo.R1offset; d < feti.R1.size(); ++d) {
        dinfo.push_back(DomainInfo((int)(domainOffset + d), koffset, feti.R1[d].nrows));
        koffset += feti.R1[d].nrows;
    }

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
//    #pragma omp parallel for
    Vector_Dense<T> _e;
    _e.size = feti.sinfo.R1size;
    _e.vals = e.vals + feti.sinfo.R1offset;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        math::blas::apply(_e, T{-1}, feti.R1[d], T{0}, feti.f[d]);
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
            for (int c = 0; c < Gt.ncols; ++c) {
                out.vals[r + Vector_Kernel<T>::offset] += Gt.vals[r * Gt.ncols + c] * in.vals[c];
            }
//            for (int c = G.rows[r]; c < G.rows[r + 1]; ++c) {
//                out.vals[r + Vector_Kernel<T>::offset] += G.vals[c] * in.vals[G.cols[c]];
//            }
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
    for (int r = 0; r < Gt.nrows; ++r) {
        for (int c = 0; c < Gt.ncols; ++c) {
            out.vals[c] += alpha * Gt.vals[r * Gt.ncols + c] * in.vals[r];
        }
//        for (int c = G.rows[r]; c < G.rows[r + 1]; ++c) {
//            out.vals[G.cols[c]] += alpha * G.vals[c] * in.vals[r];
//        }
    }
    out.synchronize();
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_applyR(const Vector_Dense<T> &in, std::vector<Vector_Dense<T> > &out)
{
    for (size_t d = 0; d < out.size(); ++d) {
        math::blas::applyT(out[d], T{1}, feti.R1[d], T{0}, in);
    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_computeDualGraph()
{
    dualGraph.resize(1);
    size_t n = 0;
    for (; n < feti.decomposition->neighbors.size() && feti.decomposition->neighbors[n] < info::mpi::rank; ++n) {
        dualGraph[0].push_back(DomainInfo(feti.decomposition->neighbors[n], 0, 0));
    }
    dualGraph[0].push_back(DomainInfo(info::mpi::rank, feti.sinfo.R1offset, feti.sinfo.R1size));
    for (; n < feti.decomposition->neighbors.size(); ++n) {
        dualGraph[0].push_back(DomainInfo(feti.decomposition->neighbors[n], 0, 0));
    }

//    std::vector<std::vector<DomainInfo> > sBuffer(feti.decomposition->neighbors.size()), rBuffer(feti.decomposition->neighbors.size());
//    for (size_t d = 0; d < dualGraph.size(); ++d) {
//        int last = -1;
//        for (size_t di = 0; di < dualGraph[d].size(); ++di) {
//            int n = feti.decomposition->noffset(dualGraph[d][di].domain);
//            if (!feti.decomposition->ismy(dualGraph[d][di].domain) && last < n) {
//                sBuffer[n].push_back(dinfo[d]);
//                last = n;
//            }
//        }
//    }
//
//    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, feti.decomposition->neighbors)) {
//        eslog::error("cannot exchange dual graph info\n");
//    }
//
//    std::unordered_map<int, DomainInfo> other;
//    for (size_t n = 0; n < rBuffer.size(); ++n) {
//        for (size_t i = 0; i < rBuffer[n].size(); ++i) {
//            other[rBuffer[n][i].domain] = rBuffer[n][i];
//        }
//    }
//
//    for (size_t d = 0; d < dinfo.size(); ++d) {
//        for (size_t i = 0; i < dualGraph[d].size(); ++i) {
//            if (!feti.decomposition->ismy(dualGraph[d][i].domain)) {
//                dualGraph[d][i] = other[dualGraph[d][i].domain];
//            }
//        }
//    }
//
//    downinfo.resize(feti.decomposition->neighbors.size());
//    std::vector<int> cOffset(dinfo.size());
//    for (size_t i = 0, offset = 0; i < feti.lambdas.cmap.size(); ) {
//        int lsize = feti.lambdas.cmap[i];
//        int domains = feti.lambdas.cmap[i + 1];
//        int last = -1;
//        for (int d = 0; d < domains; ++d) {
//            int di = feti.lambdas.cmap[i + 2 + d] - domainOffset;
//            if (feti.lambdas.cmap[i + 2 + d] < feti.decomposition->dbegin) {
//                int n = feti.decomposition->noffset(feti.lambdas.cmap[i + 2 + d]);
//                if (last < n) {
//                    for (int d2 = d; d2 < domains; ++d2) {
//                        if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d2])) {
//                            int di2 = feti.lambdas.cmap[i + 2 + d2] - domainOffset;
//                            downinfo[n][feti.lambdas.cmap[i + 2 + d2]] = dinfo[feti.lambdas.cmap[i + 2 + d2] - domainOffset];
//                            downinfo[n][feti.lambdas.cmap[i + 2 + d2]].cindices.push_back({cOffset[di2], lsize});
//                            downinfo[n][feti.lambdas.cmap[i + 2 + d2]].ncols += lsize;
//                        }
//                    }
//                    last = n;
//                }
//            }
//            if (feti.decomposition->dend <= feti.lambdas.cmap[i + 2 + d]) {
//                upinfo[feti.lambdas.cmap[i + 2 + d]] = other[feti.lambdas.cmap[i + 2 + d]];
//                upinfo[feti.lambdas.cmap[i + 2 + d]].cindices.push_back({(int)offset, lsize});
//                upinfo[feti.lambdas.cmap[i + 2 + d]].ncols += lsize;
//            }
//            if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d])) {
//                cOffset[di] += lsize;
//            }
//        }
//        i += feti.lambdas.cmap[i + 1] + 2;
//        offset += lsize;
//    }
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_setG()
{
    Gt.resize(feti.sinfo.R1size, feti.lambdas.size);
    eslog::checkpointln("FETI: SET G");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_updateG()
{
    math::set(Gt, T{0});
    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (int kr = 0; kr < dinfo[d].kernels; ++kr) {
            for (int c = 0; c < feti.B1[d].nrows; ++c) {
                for (int i = feti.B1[d].rows[c]; i < feti.B1[d].rows[c + 1]; ++i) {
                    Gt.vals[kr * Gt.ncols + feti.D2C[d][c]] -= feti.R1[d].vals[feti.R1[d].ncols * kr + feti.B1[d].cols[i]] * feti.B1[d].vals[i];
                }
            }
        }
    }
//    const FETIDecomposition *decomposition = feti.decomposition;
//
//    std::vector<std::vector<T> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
//    for (size_t n = 0; n < downinfo.size(); ++n) {
//        for (auto di = downinfo[n].cbegin(); di != downinfo[n].cend(); ++di) {
//            const NeighborDomainInfo &ndi = di->second;
//            for (int kr = 0; kr < ndi.kernels; ++kr) {
//                int roffset = G.rows[ndi.koffset - feti.sinfo.R1offset + kr];
//                for (size_t cc = 0; cc < ndi.cindices.size(); ++cc) {
//                    for (int c = 0; c < ndi.cindices[cc].count; ++c) {
//                        sBuffer[n].push_back(G.vals[roffset + ndi.cindices[cc].offset + c]);
//                    }
//                }
//            }
//        }
//    }
//
//    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
//        eslog::error("cannot exchange neighbor's G\n");
//    }
//
//    for (size_t n = 0, offset = G.rows[G.nrows]; n < rBuffer.size(); ++n) {
//        std::copy(rBuffer[n].begin(), rBuffer[n].end(), Gt.vals + offset);
//        offset += rBuffer[n].size();
//    }
    eslog::checkpointln("FETI: UPDATE G");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_setGGt()
{
    const int IDX = Indexing::CSR;

    GGtDataOffset = 0;
    GGtDataOffset = dualGraph[0].size();
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

    invGGt.resize(Gt.nrows, feti.sinfo.R1totalSize);
    eslog::checkpointln("FETI: SET GGT");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_updateGGt()
{
    const int IDX = Indexing::CSR;

    for (int kr = 0; kr < feti.sinfo.R1size; ++kr) {
        for (size_t i = 0, c = GGt.rows[feti.sinfo.R1offset + kr] - IDX; i < dualGraph[0].size(); ++i) {
            for (int kc = 0; kc < dualGraph[0][i].kernels; ++kc) {
                if (feti.sinfo.R1offset + kr <= dualGraph[0][i].koffset + kc) {
                    GGt.vals[c] = 0;
                    for (int kk = 0; kk < Gt.ncols; ++kk) {
                        GGt.vals[c] += Gt.vals[kk] * Gt.vals[kk];
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
        eye.resize(Gt.nrows, feti.sinfo.R1totalSize);
        math::set(eye, T{});
        for (int r = 0; r < Gt.nrows; ++r) {
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
        eslog::storedata(" STORE: feti/projector/{G, e, GGt, invGGt}\n");
//        math::store(G, utils::filename(utils::debugDirectory(step) + "/feti/projector", "G").c_str());
        math::store(Gt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "Gt").c_str());
//        math::store(e, utils::filename(utils::debugDirectory(step) + "/feti/projector", "e").c_str());
        math::store(GGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "GGt").c_str());
        math::store(invGGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invGGt").c_str());
    }
}

template struct HFETIOrthogonalSymmetric<double>;
template struct HFETIOrthogonalSymmetric<std::complex<double> >;

}
