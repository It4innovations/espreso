
#include "tfeti.orthogonal.symmetric.h"
#include "basis/containers/serializededata.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"

#include <vector>
#include <unordered_map>

namespace espreso {

template<typename T>
TFETIOrthogonalSymmetric<T>::TFETIOrthogonalSymmetric(FETI<T> &feti)
: Projector<T>(feti)
{
    kernel.resize(feti.R1.size());
    for (size_t d = 0, offset = 0; d < feti.R1.size(); ++d) {
        kernel[d].offset = offset;
        offset += kernel[d].size = feti.R1[d].nrows;
        Projector<T>::Kernel::total += feti.R1[d].nrows;
    }
    Projector<T>::Kernel::roffset = Projector<T>::Kernel::rsize = Projector<T>::Kernel::total;
    Projector<T>::Kernel::total = Communication::exscan(Projector<T>::Kernel::roffset);
    Vector_Kernel<T>::set(Projector<T>::Kernel::roffset, Projector<T>::Kernel::rsize, Projector<T>::Kernel::total);

    iGGtGx.resize(Projector<T>::Kernel::rsize);
    Gx.resize(Projector<T>::Kernel::total);
    e.resize(Projector<T>::Kernel::total);

    domainOffset = feti.decomposition->dbegin;
    dinfo.reserve(feti.R1.size());
    for (size_t d = 0, koffset = Projector<T>::Kernel::roffset; d < feti.R1.size(); ++d) {
        dinfo.push_back(DomainInfo((int)(domainOffset + d), koffset, feti.R1[d].nrows));
        koffset += feti.R1[d].nrows;
    }

    _computeDualGraph();
    _setG();
    _setGGt();
}

template<typename T>
TFETIOrthogonalSymmetric<T>::~TFETIOrthogonalSymmetric()
{

}

template<typename T>
void TFETIOrthogonalSymmetric<T>::update(const step::Step &step)
{
    #pragma omp parallel for
    for (size_t d = 0; d < dinfo.size(); ++d) {
        Vector_Dense<T> _e;
        _e.size = feti.R1[d].nrows;
        _e.vals = e.vals + dinfo[d].koffset;
        math::blas::apply(_e, T{-1}, feti.R1[d], T{0}, feti.f[d]);
    }
    e.synchronize();
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [e]");

    if (feti.updated.K || feti.updated.B) {
        _updateG();
        _updateGGt();
    }
    Projector<T>::_print(step);
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_computeDualGraph()
{
    dualGraph.resize(dinfo.size());
    for (size_t i = 0; i < feti.lambdas.cmap.size(); ) {
        int domains = feti.lambdas.cmap[i + 1];
        for (int d1 = 0; d1 < domains; ++d1) {
            int di1 = feti.lambdas.cmap[i + 2 + d1] - domainOffset;
            for (int d2 = 0; d2 < domains; ++d2) {
                int di2 = feti.lambdas.cmap[i + 2 + d2] - domainOffset;
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
        int last = -1;
        for (size_t di = 0; di < dualGraph[d].size(); ++di) {
            int n = feti.decomposition->noffset(dualGraph[d][di].domain);
            if (!feti.decomposition->ismy(dualGraph[d][di].domain) && last < n) {
                sBuffer[n].push_back(dinfo[d]);
                last = n;
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, feti.decomposition->neighbors)) {
        eslog::error("cannot exchange dual graph info\n");
    }

    std::unordered_map<int, DomainInfo> other;
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
    std::vector<int> cOffset(dinfo.size());
    for (size_t i = 0, offset = 0; i < feti.lambdas.cmap.size(); ) {
        int lsize = feti.lambdas.cmap[i];
        int domains = feti.lambdas.cmap[i + 1];
        int last = -1;
        for (int d = 0; d < domains; ++d) {
            int di = feti.lambdas.cmap[i + 2 + d] - domainOffset;
            if (feti.lambdas.cmap[i + 2 + d] < feti.decomposition->dbegin) {
                int n = feti.decomposition->noffset(feti.lambdas.cmap[i + 2 + d]);
                if (last < n) {
                    for (int d2 = d; d2 < domains; ++d2) {
                        if (feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d2])) {
                            int di2 = feti.lambdas.cmap[i + 2 + d2] - domainOffset;
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
                upinfo[feti.lambdas.cmap[i + 2 + d]].cindices.push_back({(int)offset, lsize});
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
void TFETIOrthogonalSymmetric<T>::_setG()
{
    // G is stored with 0-based in indexing
    int Grows = 0, Gnnz = 0;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        Grows += dinfo[d].kernels;
        Gnnz += dinfo[d].kernels * feti.D2C[d].size();
    }
    int Gtrows = Grows, Gtnnz = Gnnz;
    for (auto di = upinfo.cbegin(); di != upinfo.cend(); ++di) {
        Gtrows += di->second.kernels;
        Gtnnz += di->second.kernels * di->second.ncols;
    }

    Gt.resize(Gtrows, feti.lambdas.size, Gtnnz);
    Gt.rows[0] = 0;
    int ri = 0;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (int kr = 0; kr < dinfo[d].kernels; ++kr, ++ri) {
            Gt.rows[ri + 1] = Gt.rows[ri] + feti.B1[d].nrows;
            for (int c = 0; c < feti.B1[d].nrows; ++c) {
                Gt.cols[Gt.rows[ri] + c] = feti.D2C[d][c];
            }
        }
    }
    for (auto di = upinfo.begin(); di != upinfo.end(); ++di) {
        NeighborDomainInfo &ndi = di->second;
        ndi.koffset = ri;
        for (int kr = 0; kr < ndi.kernels; ++kr, ++ri) {
            Gt.rows[ri + 1] = Gt.rows[ri] + ndi.ncols;
            for (size_t ci = 0, c = 0; ci < ndi.cindices.size(); ++ci) {
                for (int cc = 0; cc < ndi.cindices[ci].count; ++cc, ++c) {
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
void TFETIOrthogonalSymmetric<T>::_updateG()
{
    // G is stored with 0-based in indexing
    for (size_t d = 0, r = 0; d < dinfo.size(); ++d) {
        for (int kr = 0; kr < dinfo[d].kernels; ++kr, ++r) {
            for (int c = 0; c < feti.B1[d].nrows; ++c) {
                G.vals[G.rows[r] + c] = 0;
                for (int i = feti.B1[d].rows[c]; i < feti.B1[d].rows[c + 1]; ++i) {
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
            for (int kr = 0; kr < ndi.kernels; ++kr) {
                int roffset = G.rows[ndi.koffset - Projector<T>::Kernel::roffset + kr];
                for (size_t cc = 0; cc < ndi.cindices.size(); ++cc) {
                    for (int c = 0; c < ndi.cindices[cc].count; ++c) {
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
void TFETIOrthogonalSymmetric<T>::_setGGt()
{
    const int IDX = Indexing::CSR;

    GGtDataOffset = 0;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (int kr = 0; kr < dinfo[d].kernels; ++kr) {
            for (size_t i = 0; i < dualGraph[d].size(); ++i) {
                for (int kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
                    if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
                        ++GGtDataOffset;
                    }
                }
            }
        }
    }
    GGtDataSize = GGtDataOffset;
    GGtNnz = Communication::exscan(GGtDataOffset);

    GGt.resize(Projector<T>::Kernel::total, Projector<T>::Kernel::total, GGtNnz);
    GGt.shape = Matrix_Shape::UPPER;
    GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    GGt.rows[0] = IDX;
    GGt.rows[Projector<T>::Kernel::roffset] = GGtDataOffset + IDX;
    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (int kr = 0; kr < dinfo[d].kernels; ++kr) {
            GGt.rows[dinfo[d].koffset + kr + 1] = GGt.rows[dinfo[d].koffset + kr];
            for (size_t i = 0, c = GGt.rows[dinfo[d].koffset + kr] - IDX; i < dualGraph[d].size(); ++i) {
                for (int kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
                    if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
                        GGt.cols[c++] = dualGraph[d][i].koffset + kc + IDX;
                        ++GGt.rows[dinfo[d].koffset + kr + 1];
                    }
                }
            }
        }
    }

    if (!Communication::allGatherInplace(GGt.rows, Projector<T>::Kernel::roffset + 1, G.nrows)) {
        eslog::error("cannot gather GGt rows.\n");
    }
    if (!Communication::allGatherInplace(GGt.cols, GGtDataOffset, GGtDataSize)) {
        eslog::error("cannot gather GGt cols.\n");
    }

    invGGt.resize(G.nrows, Projector<T>::Kernel::total);
    eslog::checkpointln("FETI: SET GGT");
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_updateGGt()
{
    const int IDX = Indexing::CSR;

    for (size_t d = 0; d < dinfo.size(); ++d) {
        for (int kr = 0; kr < dinfo[d].kernels; ++kr) {
            for (size_t i = 0, c = GGt.rows[dinfo[d].koffset + kr] - IDX; i < dualGraph[d].size(); ++i) {
                for (int kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
                    if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
                        GGt.vals[c] = 0;
                        int k1, k2, ke1, ke2;
                        k1  = G.rows[dinfo[d].koffset - Projector<T>::Kernel::roffset + kr];
                        ke1 = G.rows[dinfo[d].koffset - Projector<T>::Kernel::roffset + kr + 1];
                        if (dualGraph[d][i].koffset - Projector<T>::Kernel::roffset < G.nrows) {
                            k2  = G.rows[dualGraph[d][i].koffset - Projector<T>::Kernel::roffset + kc];
                            ke2 = G.rows[dualGraph[d][i].koffset - Projector<T>::Kernel::roffset + kc + 1];
                        } else {
                            k2  = Gt.rows[upinfo[dualGraph[d][i].domain].koffset + kc];
                            ke2 = Gt.rows[upinfo[dualGraph[d][i].domain].koffset + kc + 1];
                        }
                        while (k1 < ke1 && k2 < ke2) {
                            while (k1 < ke1 && Gt.cols[k1] < Gt.cols[k2]) { ++k1; };
                            while (k2 < ke2 && Gt.cols[k2] < Gt.cols[k1]) { ++k2; };
                            if (k1 < ke1 && k2 < ke2 && Gt.cols[k1] == Gt.cols[k2]) {
                                GGt.vals[c] += Gt.vals[k1++] * Gt.vals[k2++];
                            }
                        }
                        ++c;
                    }
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
        eye.resize(G.nrows, Projector<T>::Kernel::total);
        math::set(eye, T{});
        for (int r = 0; r < G.nrows; ++r) {
            eye.vals[r * Projector<T>::Kernel::total + Projector<T>::Kernel::roffset + r] = T{1};
        }
        switch (Projector<T>::GGTtype) {
        case Projector<T>::GGT_TYPE::GGT:
            GGtSolver.solve(eye, invGGt);
            break;
        case Projector<T>::GGT_TYPE::LU:
            GGtSolver.solveForward(eye, invU);
            GGtSolver.solveBackward(eye, invL);
            break;
        default: break;
        }
    }
    eslog::checkpointln("FETI: COMPUTE GGT INVERSE");
}

template struct TFETIOrthogonalSymmetric<double>;
template struct TFETIOrthogonalSymmetric<std::complex<double> >;

}
