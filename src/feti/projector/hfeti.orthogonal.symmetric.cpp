
#include "hfeti.orthogonal.symmetric.h"
#include "basis/containers/serializededata.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"

#include <vector>
#include <unordered_map>

namespace espreso {

template<typename T>
HFETIOrthogonalSymmetric<T>::HFETIOrthogonalSymmetric(FETI<T> &feti)
: Projector<T>(feti), domainOffset(0), GGtDataOffset(0), GGtDataSize(0), GGtNnz(0)
{

}

template<typename T>
HFETIOrthogonalSymmetric<T>::~HFETIOrthogonalSymmetric()
{

}

template<typename T>
void HFETIOrthogonalSymmetric<T>::set(const step::Step &step)
{

}

template<typename T>
void HFETIOrthogonalSymmetric<T>::update(const step::Step &step)
{
    kernel.resize(feti.R1.size());
    for (size_t d = 0; d < feti.R1.size(); ++d) {
        kernel[d].offset = 0;
        kernel[d].size = feti.R1[d].nrows;
        Projector<T>::Kernel::total = std::max(feti.R1[d].nrows, Projector<T>::Kernel::total);
    }
    Projector<T>::Kernel::roffset = Projector<T>::Kernel::rsize = Projector<T>::Kernel::total;
    Projector<T>::Kernel::total = Communication::exscan(Projector<T>::Kernel::roffset);
    Vector_Kernel<T>::set(Projector<T>::Kernel::roffset, Projector<T>::Kernel::rsize, Projector<T>::Kernel::total);

    iGGtGx.resize(Projector<T>::Kernel::rsize);
    Gx.resize(Projector<T>::Kernel::total);
    e.resize(Projector<T>::Kernel::total);

    cinfo.reserve(1);
    cinfo.push_back(ClusterInfo(info::mpi::rank, Projector<T>::Kernel::roffset, Projector<T>::Kernel::rsize));

    _computeDualGraph();
    _setG();
    _setGGt();

    /////

    Vector_Dense<T> _e;
    _e.size = Projector<T>::Kernel::rsize;
    _e.vals = e.vals + Projector<T>::Kernel::roffset;
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
    Projector<T>::_print(step);
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
    neighInfo.resize(feti.decomposition->neighbors.size());
    size_t n = 0;
    for (; n < feti.decomposition->neighbors.size() && feti.decomposition->neighbors[n] < info::mpi::rank; ++n) {
        dualGraph[0].push_back(ClusterInfo(feti.decomposition->neighbors[n], rBuffer[n][0].koffset, rBuffer[n][0].kernels));
        neighInfo[n] = dualGraph[0].back();
    }
    dualGraph[0].push_back(ClusterInfo(info::mpi::rank, Projector<T>::Kernel::roffset, Projector<T>::Kernel::rsize));
    for (; n < feti.decomposition->neighbors.size(); ++n) {
        dualGraph[0].push_back(ClusterInfo(feti.decomposition->neighbors[n], rBuffer[n][0].koffset, rBuffer[n][0].kernels));
        neighInfo[n] = dualGraph[0].back();
    }

    for (size_t i = 0, offset = 0; i < feti.lambdas.cmap.size(); ) {
        int lsize = feti.lambdas.cmap[i];
        int domains = feti.lambdas.cmap[i + 1];
        for (int d = 0, last = -1; d < domains; ++d) {
            int n = feti.decomposition->noffset(feti.lambdas.cmap[i + 2 + d]);
            if (!feti.decomposition->ismy(feti.lambdas.cmap[i + 2 + d]) && last < n) {
                if (feti.lambdas.cmap[i + 2 + d] < feti.decomposition->dbegin) {
                    neighInfo[n].cindices.push_back({ (int)offset, lsize });
                    neighInfo[n].ncols += lsize;
                }
                if (feti.decomposition->dend <= feti.lambdas.cmap[i + 2 + d]) {
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
    int Gtrows = Projector<T>::Kernel::rsize, Gtnnz = Projector<T>::Kernel::rsize * feti.lambdas.size;
    for (size_t n = 0; n < feti.decomposition->neighbors.size(); ++n) {
        if (info::mpi::rank < feti.decomposition->neighbors[n]) {
            Gtrows += neighInfo[n].kernels;
            Gtnnz  += neighInfo[n].kernels * neighInfo[n].ncols;
        }
    }

    Gt.resize(Gtrows, feti.lambdas.size, Gtnnz);
    Gt.rows[0] = 0;
    int ri = 0;
    for (int kr = 0; kr < Projector<T>::Kernel::rsize; ++kr, ++ri) {
        Gt.rows[ri + 1] = Gt.rows[ri] + feti.lambdas.size;
        for (int c = 0; c < feti.lambdas.size; ++c) {
            Gt.cols[Gt.rows[ri] + c] = c;
        }
    }
    for (size_t n = 0; n < feti.decomposition->neighbors.size(); ++n) {
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

    G.shallowCopy(Gt);
    G.nrows = Projector<T>::Kernel::rsize;
    G.nnz = Projector<T>::Kernel::rsize * feti.lambdas.size;
    eslog::checkpointln("FETI: SET G");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_updateG()
{
    math::set(G, T{0});
    for (size_t d = 0; d < feti.K.size(); ++d) {
        for (int kr = 0; kr < feti.R1[d].nrows; ++kr) {
            for (int c = 0; c < feti.B1[d].nrows; ++c) {
                for (int i = feti.B1[d].rows[c]; i < feti.B1[d].rows[c + 1]; ++i) {
                    G.vals[G.rows[kr] + feti.D2C[d][c]] -= feti.R1[d].vals[feti.R1[d].ncols * kr + feti.B1[d].cols[i]] * feti.B1[d].vals[i];
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
                        sBuffer[n].push_back(G.vals[G.rows[kr] + neighInfo[n].cindices[cc].offset + c]);
                    }
                }
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
        eslog::error("cannot exchange neighbor's G\n");
    }

    for (size_t n = 0, offset = Gt.rows[Projector<T>::Kernel::rsize]; n < rBuffer.size(); ++n) {
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

    GGt.resize(Projector<T>::Kernel::total, Projector<T>::Kernel::total, GGtNnz);
    GGt.shape = Matrix_Shape::UPPER;
    GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    GGt.rows[0] = IDX;
    GGt.rows[Projector<T>::Kernel::roffset] = GGtDataOffset + IDX;

    for (int kr = 0; kr < Projector<T>::Kernel::rsize; ++kr) {
        GGt.rows[Projector<T>::Kernel::roffset + kr + 1] = GGt.rows[Projector<T>::Kernel::roffset + kr];
        for (size_t i = 0, c = GGt.rows[Projector<T>::Kernel::roffset + kr] - IDX; i < dualGraph[0].size(); ++i) {
            for (int kc = 0; kc < dualGraph[0][i].kernels; ++kc) {
                if (Projector<T>::Kernel::roffset + kr <= dualGraph[0][i].koffset + kc) {
                    GGt.cols[c++] = dualGraph[0][i].koffset + kc + IDX;
                    ++GGt.rows[Projector<T>::Kernel::roffset + kr + 1];
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

    switch (Projector<T>::GGTtype) {
    case Projector<T>::GGT_TYPE::GGT:
        invGGt.resize(G.nrows, Projector<T>::Kernel::total);
        break;
    case Projector<T>::GGT_TYPE::LU:
        invU.resize(G.nrows, Projector<T>::Kernel::total);
        invL.resize(G.nrows, Projector<T>::Kernel::total);
        break;
    default: break;
    }
    eslog::checkpointln("FETI: SET GGT");
}

template<typename T>
void HFETIOrthogonalSymmetric<T>::_updateGGt()
{
    const int IDX = Indexing::CSR;
    for (int kr = 0; kr < Projector<T>::Kernel::rsize; ++kr) {
        for (size_t i = 0, c = GGt.rows[cinfo[0].koffset + kr] - IDX; i < dualGraph[0].size(); ++i) {
            for (int kc = 0; kc < dualGraph[0][i].kernels; ++kc) {
                if (Projector<T>::Kernel::roffset + kr <= dualGraph[0][i].koffset + kc) {
                    GGt.vals[c] = 0;
                    int k1, k2, ke1, ke2;
                    k1  = G.rows[kr];
                    ke1 = G.rows[kr + 1];
                    if (dualGraph[0][i].koffset - Projector<T>::Kernel::roffset < G.nrows) {
                        k2  = G.rows[kc];
                        ke2 = G.rows[kc + 1];
                    } else {
                        int koffset = neighInfo[i - 1].koffset; // (i - i) since there is local cluster in dualGraph
                        k2  = Gt.rows[koffset + kc];
                        ke2 = Gt.rows[koffset + kc + 1];
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
    } else {
        Projector<T>::GGTtype = Projector<T>::GGT_TYPE::NONE;
    }
    eslog::checkpointln("FETI: COMPUTE GGT INVERSE");
}

template struct HFETIOrthogonalSymmetric<double>;
template struct HFETIOrthogonalSymmetric<std::complex<double> >;

}
