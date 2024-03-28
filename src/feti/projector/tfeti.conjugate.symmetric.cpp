
#include "tfeti.conjugate.symmetric.h"
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
TFETIConjugateSymmetric<T>::TFETIConjugateSymmetric(FETI<T> &feti)
: Projector<T>(feti)
{
    iGFGtGx.resize(feti.sinfo.R1size);

    domainOffset = feti.decomposition->dbegin;

    dinfo.reserve(feti.R1.size());
    for (size_t d = 0, koffset = feti.sinfo.R1offset; d < feti.R1.size(); ++d) {
        dinfo.push_back(DomainInfo((esint)(domainOffset + d), koffset, feti.R1[d].nrows));
        koffset += feti.R1[d].nrows;
    }

    _computeDualGraph();
    _setG();
    _setGFGt();
}

template<typename T>
TFETIConjugateSymmetric<T>::~TFETIConjugateSymmetric()
{

}

template<typename T>
void TFETIConjugateSymmetric<T>::info()
{
    esint nnz = 2 * (GFGt.nnz - GFGt.nrows) + GFGt.nrows;

    eslog::info(" = CONJUGATE PROJECTOR PROPERTIES                                                            = \n");
    eslog::info(" =   GFGT ROWS                                                                     %9d = \n", GFGt.nrows);
    eslog::info(" =   GFGT NNZ                                                                      %9d = \n", nnz);
//    eslog::info(" =   GFGt FACTORS NNZ                                                               %9d = \n", GFGtSolver.nnzL);
    if (feti.configuration.exhaustive_info) {
        // PPt = eye
    }
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template<typename T>
void TFETIConjugateSymmetric<T>::update(const step::Step &step)
{
    #pragma omp parallel for
    for (size_t d = 0; d < dinfo.size(); ++d) {
        if (feti.updated.K) {
            math::orthonormalize(feti.R1[d]);
        }
        Vector_Dense<T> _e;
        _e.size = feti.R1[d].nrows;
        _e.vals = e.vals + dinfo[d].koffset;
        math::blas::apply(_e, T{1}, feti.R1[d], T{0}, feti.f[d]);
    }
    e.synchronize();
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [e]");

    if (feti.updated.K || feti.updated.B) {
        _updateG();
        _updateGFGt();
    }
    _print(step);
}

template<typename T>
void TFETIConjugateSymmetric<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    if (GFGt.nrows) {
        x.copyToWithoutHalo(y);
        _applyG(x, Gx);
        _applyInvGFGt(Gx, iGFGtGx);
        _applyGt(iGFGtGx, T{-1}, y);
    } else {
        math::copy(y, x);
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::apply_e(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
    if (GFGt.nrows) {
        math::set(y, T{0});
        _applyInvGFGt(x, iGFGtGx);
        _applyGt(iGFGtGx, T{-1}, y);
    } else {
        math::set(y, T{0});
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::apply_R(const Vector_Kernel<T> &x,  std::vector<Vector_Dense<T> > &y)
{
    if (GFGt.nrows) {
        _applyR(x, y);
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < y.size(); ++d) {
            math::set(y[d], T{0});
        }
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::apply_Ra(const Vector_Dual<T> &x,  std::vector<Vector_Dense<T> > &y)
{
    if (GFGt.nrows) {
        _applyG(x, Gx);
        _applyInvGFGt(Gx, iGFGtGx);
        _applyR(iGFGtGx, y);
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < y.size(); ++d) {
            math::set(y[d], T{0});
        }
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::apply_invU(const Vector_Kernel<T> &x, Vector_Kernel<T> &y)
{

}

template<typename T>
void TFETIConjugateSymmetric<T>::apply_invL(const Vector_Kernel<T> &x, Vector_Kernel<T> &y)
{

}


template<typename T>
void TFETIConjugateSymmetric<T>::apply_GtinvU(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{

}

template<typename T>
void TFETIConjugateSymmetric<T>::apply_invLG(const Vector_Dual<T> &x, Vector_Kernel<T> &y)
{

}

template<typename T>
void TFETIConjugateSymmetric<T>::_applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out)
{
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t r = Vector_Kernel<T>::distribution[t]; r < Vector_Kernel<T>::distribution[t + 1]; ++r) {
            out.vals[r + Vector_Kernel<T>::offset] = T{0};
            for (esint c = G.rows[r]; c < G.rows[r + 1]; ++c) {
                out.vals[r + Vector_Kernel<T>::offset] += G.vals[c] * in.vals[G.cols[c]];
            }
        }
    }
    out.synchronize();
}

template<typename T>
void TFETIConjugateSymmetric<T>::_applyInvGFGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        Matrix_Dense<T> a;
        Vector_Dense<T> y;
        a.ncols = invGFGt.ncols;
        a.nrows = y.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];

        a.vals = invGFGt.vals + invGFGt.ncols * Vector_Kernel<T>::distribution[t];
        y.vals = out.vals + Vector_Kernel<T>::distribution[t];

        math::blas::apply(y, T{1}, a, T{0}, in);
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::_applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out)
{
    for (esint r = 0; r < G.nrows; ++r) {
        for (esint c = G.rows[r]; c < G.rows[r + 1]; ++c) {
            out.vals[G.cols[c]] += alpha * G.vals[c] * in.vals[r];
        }
    }
    out.synchronize();
}

template<typename T>
void TFETIConjugateSymmetric<T>::_applyR(const Vector_Dense<T> &in, std::vector<Vector_Dense<T> > &out)
{
    #pragma omp parallel for
    for (size_t d = 0; d < out.size(); ++d) {
        Vector_Dense<T> y;
        y.size = dinfo[d].kernels;
        y.vals = in.vals + dinfo[d].koffset - feti.sinfo.R1offset;

        math::blas::applyT(out[d], T{1}, feti.R1[d], T{0}, y);
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::_computeDualGraph()
{
    dualGraph.resize(dinfo.size());
    for (size_t i = 0; i < feti.lambdas.cmap.size(); ) {
        esint domains = feti.lambdas.cmap[i + 1];
        for (esint d1 = 0; d1 < domains; ++d1) {
            esint di1 = feti.lambdas.cmap[i + 2 + d1] - domainOffset;
            for (esint d2 = 0; d2 < domains; ++d2) {
                esint di2 = feti.lambdas.cmap[i + 2 + d2] - domainOffset;
                if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d1]) && feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d2])) {
                    dualGraph[di1].push_back(dinfo[di2]);
                } else {
                    if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d1]) && !feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d2])) {
                        dualGraph[di1].push_back(DomainInfo(feti.lambdas.cmap[i + 2 + d2], 0, 0));
                    }
                    if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d2]) && !feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d1])) {
                        dualGraph[di2].push_back(DomainInfo(feti.lambdas.cmap[i + 2 + d1], 0, 0));
                    }
                }
            }
        }
        i += feti.lambdas.cmap[i + 1] + 2;
    }
    for (size_t d = 0; d < dinfo.size(); ++d) {
        utils::sortAndRemoveDuplicates(dualGraph[d]);
    }

    std::vector<std::vector<DomainInfo> > sBuffer(feti.decomposition->neighbors.size()), rBuffer(feti.decomposition->neighbors.size());
    for (size_t d = 0; d < dualGraph.size(); ++d) {
        esint last = -1;
        for (size_t di = 0; di < dualGraph[d].size(); ++di) {
            esint n = feti.decomposition->noffset(dualGraph[d][di].domain);
            if (!feti.decomposition->ismy(dualGraph[d][di].domain) && last < n) {
                sBuffer[n].push_back(dinfo[d]);
                last = n;
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, feti.decomposition->neighbors)) {
        eslog::error("cannot exchange dual graph info\n");
    }

    std::unordered_map<esint, DomainInfo> other;
    for (size_t n = 0; n < rBuffer.size(); ++n) {
        for (size_t i = 0; i < rBuffer[n].size(); ++i) {
            other[rBuffer[n][i].domain] = rBuffer[n][i];
        }
    }

    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (size_t i = 0; i < dualGraph[d].size(); ++i) {
            if (!feti.decomposition->ismy(dualGraph[d][i].domain)) {
                dualGraph[d][i] = other[dualGraph[d][i].domain];
            }
        }
    }

    downinfo.resize(feti.decomposition->neighbors.size());
    std::vector<esint> cOffset(dinfo.size());
    for (size_t i = 0, offset = 0; i < feti.lambdas.cmap.size(); ) {
        esint lsize = feti.lambdas.cmap[i];
        esint domains = feti.lambdas.cmap[i + 1];
        esint last = -1;
        for (esint d = 0; d < domains; ++d) {
            esint di = feti.lambdas.cmap[i + 2 + d] - domainOffset;
            if (feti.lambdas.cmap[i + 2 + d] < feti.decomposition->dbegin) {
                esint n = feti.decomposition->noffset(feti.lambdas.cmap[i + 2 + d]);
                if (last < n) {
                    for (esint d2 = d; d2 < domains; ++d2) {
                        if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d2])) {
                            esint di2 = feti.lambdas.cmap[i + 2 + d2] - domainOffset;
                            downinfo[n][feti.lambdas.cmap[i + 2 + d2]] = dinfo[feti.lambdas.cmap[i + 2 + d2] - domainOffset];
                            downinfo[n][feti.lambdas.cmap[i + 2 + d2]].cindices.push_back({cOffset[di2], lsize});
                            downinfo[n][feti.lambdas.cmap[i + 2 + d2]].ncols += lsize;
                        }
                    }
                    last = n;
                }
            }
            if (feti.decomposition->dend <= feti.lambdas.cmap[i + 2 + d]) {
                upinfo[feti.lambdas.cmap[i + 2 + d]] = other[feti.lambdas.cmap[i + 2 + d]];
                upinfo[feti.lambdas.cmap[i + 2 + d]].cindices.push_back({(esint)offset, lsize});
                upinfo[feti.lambdas.cmap[i + 2 + d]].ncols += lsize;
            }
            if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d])) {
                cOffset[di] += lsize;
            }
        }
        i += feti.lambdas.cmap[i + 1] + 2;
        offset += lsize;
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::_setG()
{
    // G is stored with 0-based in indexing
    esint Grows = 0, Gnnz = 0;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        Grows += dinfo[d].kernels;
        Gnnz += dinfo[d].kernels * feti.D2C[d].size();
    }
    esint Gtrows = Grows, Gtnnz = Gnnz;
    for (auto di = upinfo.cbegin(); di != upinfo.cend(); ++di) {
        Gtrows += di->second.kernels;
        Gtnnz += di->second.kernels * di->second.ncols;
    }

    Gt.resize(Gtrows, feti.lambdas.size, Gtnnz);
    Gt.rows[0] = 0;
    size_t ri = 0;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (esint kr = 0; kr < dinfo[d].kernels; ++kr, ++ri) {
            Gt.rows[ri + 1] = Gt.rows[ri] + feti.B1[d].nrows;
            for (esint c = 0; c < feti.B1[d].nrows; ++c) {
                Gt.cols[Gt.rows[ri] + c] = feti.D2C[d][c];
            }
        }
    }
    for (auto di = upinfo.begin(); di != upinfo.end(); ++di) {
        NeighborDomainInfo &ndi = di->second;
        ndi.koffset = ri;
        for (esint kr = 0; kr < ndi.kernels; ++kr, ++ri) {
            Gt.rows[ri + 1] = Gt.rows[ri] + ndi.ncols;
            for (size_t ci = 0, c = 0; ci < ndi.cindices.size(); ++ci) {
                for (esint cc = 0; cc < ndi.cindices[ci].count; ++cc, ++c) {
                    Gt.cols[Gt.rows[ri] + c] = ndi.cindices[ci].offset + cc;
                }
            }
        }
    }


    G.shallowCopy(Gt);
    G.nrows = Grows;
    G.nnz = Gnnz;
    eslog::checkpointln("FETI: SET G");
}

template<typename T>
void TFETIConjugateSymmetric<T>::_updateG()
{
    // G is stored with 0-based in indexing
    for (size_t d = 0, r = 0; d < dinfo.size(); ++d) {
        for (esint kr = 0; kr < dinfo[d].kernels; ++kr, ++r) {
            for (esint c = 0; c < feti.B1[d].nrows; ++c) {
                G.vals[G.rows[r] + c] = 0;
                for (esint i = feti.B1[d].rows[c]; i < feti.B1[d].rows[c + 1]; ++i) {
                    G.vals[G.rows[r] + c] -= feti.R1[d].vals[feti.R1[d].ncols * kr + feti.B1[d].cols[i]] * feti.B1[d].vals[i];
                }
            }
        }
    }

    const FETIDecomposition *decomposition = feti.decomposition;

    std::vector<std::vector<T> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
    for (size_t n = 0; n < downinfo.size(); ++n) {
        for (auto di = downinfo[n].cbegin(); di != downinfo[n].cend(); ++di) {
            const NeighborDomainInfo &ndi = di->second;
            for (esint kr = 0; kr < ndi.kernels; ++kr) {
                esint roffset = G.rows[ndi.koffset - feti.sinfo.R1offset + kr];
                for (size_t cc = 0; cc < ndi.cindices.size(); ++cc) {
                    for (esint c = 0; c < ndi.cindices[cc].count; ++c) {
                        sBuffer[n].push_back(G.vals[roffset + ndi.cindices[cc].offset + c]);
                    }
                }
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
        eslog::error("cannot exchange neighbor's G\n");
    }

    for (size_t n = 0, offset = G.rows[G.nrows]; n < rBuffer.size(); ++n) {
        std::copy(rBuffer[n].begin(), rBuffer[n].end(), Gt.vals + offset);
        offset += rBuffer[n].size();
    }

    eslog::checkpointln("FETI: UPDATE G");
}

template<typename T>
void TFETIConjugateSymmetric<T>::_setGFGt()
{
    const int IDX = Indexing::CSR;

    GFGtDataOffset = 0;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (esint kr = 0; kr < dinfo[d].kernels; ++kr) {
            for (size_t i = 0; i < dualGraph[d].size(); ++i) {
                for (esint kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
                    if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
                        ++GFGtDataOffset;
                    }
                }
            }
        }
    }
    GFGtDataSize = GFGtDataOffset;
    GFGtNnz = Communication::exscan(GFGtDataOffset);

    GFGt.resize(feti.sinfo.R1totalSize, feti.sinfo.R1totalSize, GFGtNnz);
    GFGt.shape = Matrix_Shape::UPPER;
    GFGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    GFGt.rows[0] = IDX;
    GFGt.rows[feti.sinfo.R1offset] = GFGtDataOffset + IDX;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (esint kr = 0; kr < dinfo[d].kernels; ++kr) {
            GFGt.rows[dinfo[d].koffset + kr + 1] = GFGt.rows[dinfo[d].koffset + kr];
            for (size_t i = 0, c = GFGt.rows[dinfo[d].koffset + kr] - IDX; i < dualGraph[d].size(); ++i) {
                for (esint kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
                    if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
                        GFGt.cols[c++] = dualGraph[d][i].koffset + kc + IDX;
                        ++GFGt.rows[dinfo[d].koffset + kr + 1];
                    }
                }
            }
        }
    }

    if (!Communication::allGatherInplace(GFGt.rows, feti.sinfo.R1offset + 1, G.nrows)) {
        eslog::error("cannot gather GFGt rows.\n");
    }
    if (!Communication::allGatherInplace(GFGt.cols, GFGtDataOffset, GFGtDataSize)) {
        eslog::error("cannot gather GFGt cols.\n");
    }

    invGFGt.resize(G.nrows, feti.sinfo.R1totalSize);
    eslog::checkpointln("FETI: SET GFGt");
}

template<typename T>
void TFETIConjugateSymmetric<T>::_updateGFGt()
{
    const int IDX = Indexing::CSR;

    Matrix_Dense<T> dGt;

    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (esint kr = 0; kr < dinfo[d].kernels; ++kr) {
            for (size_t i = 0, c = GFGt.rows[dinfo[d].koffset + kr] - IDX; i < dualGraph[d].size(); ++i) {
                for (esint kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
                    if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
                        GFGt.vals[c] = 0;
                        esint k1, k2, ke1, ke2;
                        k1  = G.rows[dinfo[d].koffset - feti.sinfo.R1offset + kr];
                        ke1 = G.rows[dinfo[d].koffset - feti.sinfo.R1offset + kr + 1];
                        if (dualGraph[d][i].koffset - feti.sinfo.R1offset < G.nrows) {
                            k2  = G.rows[dualGraph[d][i].koffset - feti.sinfo.R1offset + kc];
                            ke2 = G.rows[dualGraph[d][i].koffset - feti.sinfo.R1offset + kc + 1];
                        } else {
                            k2  = Gt.rows[upinfo[dualGraph[d][i].domain].koffset + kc];
                            ke2 = Gt.rows[upinfo[dualGraph[d][i].domain].koffset + kc + 1];
                        }
                        while (k1 < ke1 && k2 < ke2) {
                            while (k1 < ke1 && Gt.cols[k1] < Gt.cols[k2]) { ++k1; };
                            while (k2 < ke2 && Gt.cols[k2] < Gt.cols[k1]) { ++k2; };
                            if (k1 < ke1 && k2 < ke2 && Gt.cols[k1] == Gt.cols[k2]) {
                                GFGt.vals[c] += Gt.vals[k1++] * Gt.vals[k2++];
                            }
                        }
                        ++c;
                    }
                }
            }
        }
    }


    if (!Communication::allGatherInplace(GFGt.vals, GFGtDataOffset, GFGtDataSize)) {
        eslog::error("cannot gather GFGt vals.\n");
    }
    eslog::checkpointln("FETI: GATHER GFGt VALUES");

    if (GFGt.nrows) {
        DirectSparseSolver<T> GFGtSolver;
        GFGtSolver.commit(GFGt);
        GFGtSolver.symbolicFactorization();
        GFGtSolver.numericalFactorization();
        eslog::checkpointln("FETI: GFGt FACTORIZATION");

        Matrix_Dense<T> eye;
        eye.resize(G.nrows, feti.sinfo.R1totalSize);
        math::set(eye, T{});
        for (esint r = 0; r < G.nrows; ++r) {
            eye.vals[r * feti.sinfo.R1totalSize + feti.sinfo.R1offset + r] = T{1};
        }
        GFGtSolver.solve(eye, invGFGt);
    }
    eslog::checkpointln("FETI: COMPUTE GFGt INVERSE");
}

template<typename T>
void TFETIConjugateSymmetric<T>::_print(const step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/projector/{R_orth, G, e, GFGt, invGFGt}\n");
        for (size_t di = 0; di < feti.R1.size(); ++di) {
            math::store(feti.R1[di], utils::filename(utils::debugDirectory(step) + "/feti/projector", (std::string("R_orth") + std::to_string(di)).c_str()).c_str());
        }
        math::store(G, utils::filename(utils::debugDirectory(step) + "/feti/projector", "G").c_str());
        math::store(Gt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "Gt").c_str());
        // math::store(e, utils::filename(utils::debugDirectory(step) + "/feti/projector", "e").c_str());
        // math::store(GFGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "GFGt").c_str());
        // math::store(invGFGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invGFGt").c_str());
    }
}

template struct TFETIConjugateSymmetric<double>;
template struct TFETIConjugateSymmetric<std::complex<double> >;

}