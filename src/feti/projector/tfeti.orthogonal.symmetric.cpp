
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
: Projector<T>(feti), GGtDataOffset(0), GGtDataSize(0)
{

}

template<typename T>
TFETIOrthogonalSymmetric<T>::~TFETIOrthogonalSymmetric()
{

}

template<typename T>
void TFETIOrthogonalSymmetric<T>::set(const step::Step &step)
{

}

template<typename T>
void TFETIOrthogonalSymmetric<T>::update(const step::Step &step)
{
    int reset = kernel.size() != feti.R1.size() || Gt.ncols != feti.lambdas.size;
    for (size_t i = 0; !reset && i < feti.R1.size(); ++i) {
        reset |= kernel[i].size != feti.R1[i].nrows;
    }
    Communication::allReduce(&reset, nullptr, 1, MPITools::getType(reset).mpitype, MPI_MAX);

    if (reset) { // only different number of kernels
        Projector<T>::reset();
        GGtDataOffset = GGtDataSize = 0;
        kernel.clear();
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

        for (size_t d = 0; d < feti.R1.size(); ++d) {
            dual.pushVertex(feti.decomposition->dbegin + d, feti.R1[d].nrows);
        }
        dual.initVertices();
        dual.setFromDomains(feti.decomposition, feti.lambdas.cmap);

        _setG();
        _setGGt();
    }

    #pragma omp parallel for
    for (size_t d = 0; d < feti.R1.size(); ++d) {
        Vector_Dense<T> _e;
        _e.size = feti.R1[d].nrows;
        _e.vals = e.vals + dual.vertices.at(feti.decomposition->dbegin + d).kernel.goffset;
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
void TFETIOrthogonalSymmetric<T>::_setG()
{
    // G is stored with 0-based in indexing
    auto vbegin = dual.vertices.find(feti.decomposition->dbegin);

    int Grows = 0, Gnnz = 0;
    int Gtrows = 0, Gtnnz = 0;
    for (auto v = vbegin; v != dual.vertices.cend(); ++v) {
        v->second.kernel.loffset = Gtrows;
        Gtrows += v->second.kernel.size;
        Gtnnz += v->second.kernel.size * v->second.lambdas.total;
        if (v->first < feti.decomposition->dend) {
            Grows += v->second.kernel.size;
            Gnnz += v->second.kernel.size * v->second.lambdas.total;
        }
    }

    Gt.resize(Gtrows, feti.lambdas.size, Gtnnz);
    Gt.rows[0] = 0;
    int ri = 0;
    for (auto v = vbegin; v != dual.vertices.cend(); ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr, ++ri) {
            Gt.rows[ri + 1] = Gt.rows[ri] + v->second.lambdas.total;
            for (size_t ci = 0, c = Gt.rows[ri]; ci < v->second.lambdas.indices.size(); ++ci) {
                for (int cc = 0; cc < v->second.lambdas.indices[ci].size; ++c, ++cc) {
                    Gt.cols[c] = v->second.lambdas.indices[ci].offset + cc;
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
    auto vbegin = dual.vertices.find(feti.decomposition->dbegin);
    auto vend   = dual.vertices.lower_bound(feti.decomposition->dend);

    int ri = 0, d = 0;
    for (auto v = vbegin; v != vend; ++v, ++d) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr, ++ri) {
            for (int c = 0; c < feti.B1[d].nrows; ++c) {
                G.vals[G.rows[ri] + c] = 0;
                for (int i = feti.B1[d].rows[c]; i < feti.B1[d].rows[c + 1]; ++i) {
                    G.vals[G.rows[ri] + c] -= feti.R1[d].vals[feti.R1[d].ncols * kr + feti.B1[d].cols[i]] * feti.B1[d].vals[i];
                }
            }
        }
    }

    const DecompositionFETI *decomposition = feti.decomposition;
    std::vector<std::vector<T> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());

    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr, ++ri) {
            auto ci = Gt.rows[v->second.kernel.goffset + kr - Projector<T>::Kernel::roffset];
            for (size_t i = 0; i < v->second.lambdas.indices.size(); ++i) {
                for (size_t n = 0; n < v->second.lambdas.indices[i].lower_neighs.size(); ++n) {
                    for (int c = 0; c < v->second.lambdas.indices[i].size; ++c) {
                        sBuffer[v->second.lambdas.indices[i].lower_neighs[n]].push_back(Gt.vals[ci + c]);
                    }
                }
                ci += v->second.lambdas.indices[i].size;
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
    auto vbegin = dual.vertices.find(feti.decomposition->dbegin);
    auto vend   = dual.vertices.lower_bound(feti.decomposition->dend);

    GGtDataOffset = 0;
    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr) {
            for (auto e = dual.edges[v->first].cbegin(); e != dual.edges[v->first].cend(); ++e) {
                for (int kc = 0; kc < dual.vertices[*e].kernel.size; ++kc) {
                    if (v->second.kernel.goffset + kr <= dual.vertices[*e].kernel.goffset + kc) {
                        ++GGtDataOffset;
                    }
                }
            }
        }
    }
    GGtDataSize = GGtDataOffset;
    size_t GGtNnz = Communication::exscan(GGtDataOffset);

    GGt.resize(Projector<T>::Kernel::total, Projector<T>::Kernel::total, GGtNnz);
    GGt.shape = Matrix_Shape::UPPER;
    GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    GGt.rows[0] = IDX;
    GGt.rows[Projector<T>::Kernel::roffset] = GGtDataOffset + IDX;
    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr) {
            GGt.rows[v->second.kernel.goffset + kr + 1] = GGt.rows[v->second.kernel.goffset + kr];
            int c = GGt.rows[v->second.kernel.goffset + kr] - IDX;
            for (auto e = dual.edges[v->first].cbegin(); e != dual.edges[v->first].cend(); ++e) {
                for (int kc = 0; kc < dual.vertices[*e].kernel.size; ++kc) {
                    if (v->second.kernel.goffset + kr <= dual.vertices[*e].kernel.goffset + kc) {
                        GGt.cols[c++] = dual.vertices[*e].kernel.goffset + kc + IDX;
                        ++GGt.rows[v->second.kernel.goffset + kr + 1];
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
void TFETIOrthogonalSymmetric<T>::_updateGGt()
{
    const int IDX = Indexing::CSR;
    auto vbegin = dual.vertices.find(feti.decomposition->dbegin);
    auto vend   = dual.vertices.lower_bound(feti.decomposition->dend);

    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr) {
            int c = GGt.rows[v->second.kernel.goffset + kr] - IDX;
            for (auto e = dual.edges[v->first].cbegin(); e != dual.edges[v->first].cend(); ++e) {
                for (int kc = 0; kc < dual.vertices[*e].kernel.size; ++kc) {
                    if (v->second.kernel.goffset + kr <= dual.vertices[*e].kernel.goffset + kc) {
                        int r1 = v->second.kernel.loffset + kr;
                        int r2 = dual.vertices[*e].kernel.loffset + kc;
                        GGt.vals[c] = 0;
                        int k1 = Gt.rows[r1], ke1 = Gt.rows[r1 + 1];
                        int k2 = Gt.rows[r2], ke2 = Gt.rows[r2 + 1];
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
    } else {
        Projector<T>::GGTtype = Projector<T>::GGT_TYPE::NONE;
    }
    eslog::checkpointln("FETI: COMPUTE GGT INVERSE");
}

template struct TFETIOrthogonalSymmetric<double>;
template struct TFETIOrthogonalSymmetric<std::complex<double> >;

}
