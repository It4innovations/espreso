
#include "tfeti.conjugate.symmetric.h"
#include "basis/containers/serializededata.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "feti/dualoperator/dualoperator.h"
#include "feti/common/matrix_dual_sparse.h"
#include "math/math.h"

#include <vector>
#include <numeric>

namespace espreso {

template<typename T>
TFETIConjugateSymmetric<T>::TFETIConjugateSymmetric(FETI<T> &feti)
: Projector<T>(feti), GGtDataOffset(0), GGtDataSize(0)
{

}

template<typename T>
TFETIConjugateSymmetric<T>::~TFETIConjugateSymmetric()
{

}

template<typename T>
void TFETIConjugateSymmetric<T>::set(const step::Step &step)
{

}

template<typename T>
void TFETIConjugateSymmetric<T>::update(const step::Step &step)
{
    int reset = kernel.size() != feti.KR1.size() || G.ncols != feti.lambdas.size;
    for (size_t i = 0; !reset && i < feti.KR1.size(); ++i) {
        reset |= kernel[i].size != feti.KR1[i].nrows;
    }
    Communication::allReduce(&reset, nullptr, 1, MPITools::getType(reset).mpitype, MPI_MAX);

    if (reset) { // only different number of kernels
        Projector<T>::reset();
        GGtDataOffset = GGtDataSize = 0;
        kernel.clear();
        kernel.resize(feti.KR1.size());
        for (size_t d = 0, offset = 0; d < feti.KR1.size(); ++d) {
            kernel[d].offset = offset;
            offset += kernel[d].size = feti.KR1[d].nrows;
            Projector<T>::Kernel::total += feti.KR1[d].nrows;
        }
        Projector<T>::Kernel::roffset = Projector<T>::Kernel::rsize = Projector<T>::Kernel::total;
        Projector<T>::Kernel::total = Communication::exscan(Projector<T>::Kernel::roffset);
        Vector_Kernel<T>::set(Projector<T>::Kernel::roffset, Projector<T>::Kernel::rsize, Projector<T>::Kernel::total);

        iGGtGx.resize(Projector<T>::Kernel::rsize);
        Gx.resize(Projector<T>::Kernel::total);

        dual.clear();
        for (size_t d = 0; d < feti.KR1.size(); ++d) {
            dual.pushVertex(d, feti.decomposition->dbegin + d, feti.KR1[d].nrows);
        }
        dual.set(feti.decomposition, feti.lambdas.cmap);
        auto singleHop = dual.domains.edges;
        dual.spread(feti.decomposition);

        nonzeros.clear();
        distributed.clear();
        std::vector<std::vector<int> > domainToRHS;
        domainToRHS.resize(feti.K.size());
        auto vbegin = dual.domains.vertices.find(feti.decomposition->dbegin);
        auto vend   = dual.domains.vertices.lower_bound(feti.decomposition->dend);
        for (auto v = vbegin; v != vend; ++v) {
            for (int k = 0; k < v->second.kernel.size; ++k) {
                nonzeros.push_back(v->second.kernel.goffset + k);
                domainToRHS[v->second.offset].push_back(v->second.kernel.goffset + k);
            }
            bool dist = false;
            for (auto e = dual.domains.edges[v->first].cbegin(); e != dual.domains.edges[v->first].cend(); ++e) {
                if (dual.domains.vertices.at(*e).rank != info::mpi::rank) {
                    for (int k = 0; k < dual.domains.vertices.at(*e).kernel.size; ++k) {
                        nonzeros.push_back(dual.domains.vertices.at(*e).kernel.goffset + k);
                        distributed.push_back(dual.domains.vertices.at(*e).kernel.goffset + k);
                    }
                    dist = true;
                }
            }
            if (dist) {
                for (int k = 0; k < v->second.kernel.size; ++k) {
                    distributed.push_back(v->second.kernel.goffset + k);
                }
            }
            for (auto e = singleHop[v->first].cbegin(); e != singleHop[v->first].cend(); ++e) {
                for (int k = 0; k < dual.domains.vertices.at(*e).kernel.size; ++k) {
                    domainToRHS[v->second.offset].push_back(dual.domains.vertices.at(*e).kernel.goffset + k);
                }
            }
        }
        utils::sortAndRemoveDuplicates(nonzeros);
        utils::sortAndRemoveDuplicates(distributed);
        filter.clear();
        filter.resize(nonzeros.size());
        for (size_t i = 0; i < domainToRHS.size(); ++i) {
            utils::sortAndRemoveDuplicates(domainToRHS[i]);
            for (size_t j = 0; j < domainToRHS[i].size(); ++j) {
                int index = std::lower_bound(nonzeros.begin(), nonzeros.end(), domainToRHS[i][j]) - nonzeros.begin();
                filter[index].push_back(i);
            }
        }

        _setG();
        _setGGt();
    }

    if (feti.updated.K || feti.updated.B) {
        _updateG();
        _updateGGt();
    }
    Projector<T>::_print(step);
}

template<typename T>
void TFETIConjugateSymmetric<T>::orthonormalizeKernels(const step::Step &step)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.KR1.size(); ++d) {
        math::orthonormalize(feti.KR1[d]);
    }
}

template<typename T>
void TFETIConjugateSymmetric<T>::_setG()
{
    // G is stored with 0-based in indexing
    auto vbegin = dual.domains.vertices.find(feti.decomposition->dbegin);
    auto vend   = dual.domains.vertices.lower_bound(feti.decomposition->dend);

    int rows = 0, nnz = 0;
    for (auto v = vbegin; v != vend; ++v) {
        rows += v->second.kernel.size;
        nnz += v->second.kernel.size * v->second.lambdas.total;
    }

    G.resize(rows, feti.lambdas.size, nnz);
    G.rows[0] = 0;
    int ri = 0;
    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr, ++ri) {
            G.rows[ri + 1] = G.rows[ri] + v->second.lambdas.total;
            for (size_t ci = 0, c = G.rows[ri]; ci < v->second.lambdas.indices.size(); ++ci) {
                for (int cc = 0; cc < v->second.lambdas.indices[ci].size; ++c, ++cc) {
                    G.cols[c] = v->second.lambdas.indices[ci].offset + cc;
                }
            }
        }
    }
    Gt.shallowCopy(G);
    eslog::checkpointln("FETI: SET G");
}

template<typename T>
void TFETIConjugateSymmetric<T>::_updateG()
{
    auto vbegin = dual.domains.vertices.find(feti.decomposition->dbegin);
    auto vend   = dual.domains.vertices.lower_bound(feti.decomposition->dend);

    int ri = 0, d = 0;
    for (auto v = vbegin; v != vend; ++v, ++d) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr, ++ri) {
            for (int c = 0; c < feti.B1[d].nrows; ++c) {
                G.vals[G.rows[ri] + c] = 0;
                for (int i = feti.B1[d].rows[c]; i < feti.B1[d].rows[c + 1]; ++i) {
                    G.vals[G.rows[ri] + c] -= feti.KR1[d].vals[feti.KR1[d].ncols * kr + feti.B1[d].cols[i]] * feti.B1[d].vals[i];
                }
            }
        }
    }

    eslog::checkpointln("FETI: UPDATE G");
}

template<typename T>
void TFETIConjugateSymmetric<T>::_setGGt()
{
    auto vbegin = dual.domains.vertices.find(feti.decomposition->dbegin);
    auto vend   = dual.domains.vertices.lower_bound(feti.decomposition->dend);

    GGtDataOffset = 0;
    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr) {
            for (auto e = dual.domains.edges[v->first].cbegin(); e != dual.domains.edges[v->first].cend(); ++e) {
                for (int kc = 0; kc < dual.domains.vertices[*e].kernel.size; ++kc) {
                    if (v->second.kernel.goffset + kr <= dual.domains.vertices[*e].kernel.goffset + kc) {
                        ++GGtDataOffset;
                    }
                }
            }
        }
    }
    GGtDataSize = GGtDataOffset;
    size_t GGtNnz = Communication::exscan(GGtDataOffset);

    const int IDX = Indexing::CSR;
    GGt.resize(Projector<T>::Kernel::total, Projector<T>::Kernel::total, GGtNnz);
    GGt.shape = Matrix_Shape::UPPER;
    GGt.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE; // set it according to type of K?
    GGt.rows[0] = IDX;
    GGt.rows[Projector<T>::Kernel::roffset] = GGtDataOffset + IDX;
    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr) {
            GGt.rows[v->second.kernel.goffset + kr + 1] = GGt.rows[v->second.kernel.goffset + kr];
            int c = GGt.rows[v->second.kernel.goffset + kr] - IDX;
            for (auto e = dual.domains.edges[v->first].cbegin(); e != dual.domains.edges[v->first].cend(); ++e) {
                for (int kc = 0; kc < dual.domains.vertices[*e].kernel.size; ++kc) {
                    if (v->second.kernel.goffset + kr <= dual.domains.vertices[*e].kernel.goffset + kc) {
                        GGt.cols[c++] = dual.domains.vertices[*e].kernel.goffset + kc + IDX;
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
void TFETIConjugateSymmetric<T>::_updateGGt()
{
    auto vbegin = dual.domains.vertices.find(feti.decomposition->dbegin);
    auto vend   = dual.domains.vertices.lower_bound(feti.decomposition->dend);

    Matrix_Dual_Sparse<T> dG(nonzeros, distributed), dFG(nonzeros, distributed);
    Vector_Dense<T> vG; vG.size = dG.ncols;
    for (size_t i = 0; i < nonzeros.size(); ++i) {
        vG.vals = dG.vals + dG.ncols * i;
        math::set(vG, T{0});
        if (Projector<T>::Kernel::roffset <= nonzeros[i] && nonzeros[i] < Projector<T>::Kernel::roffset + Projector<T>::Kernel::rsize) {
            int r = nonzeros[i] - Projector<T>::Kernel::roffset;
            for (int c = G.rows[r]; c < G.rows[r + 1]; ++c) {
                vG.vals[G.cols[c]] = G.vals[c];
            }
        }
    }
    dG.synchronize();
    feti.dualOperator->apply(dG, dFG, filter);

//    math::store((Matrix_Dense<T>)dG , ("dG" + std::to_string(info::mpi::rank)).c_str());
//    math::store((Matrix_Dense<T>)dFG, ("dFG" + std::to_string(info::mpi::rank)).c_str());

    int r1 = 0;
    for (auto v = vbegin; v != vend; ++v) {
        for (int kr = 0; kr < v->second.kernel.size; ++kr, ++r1) {
            int c = GGt.rows[v->second.kernel.goffset + kr] - Indexing::CSR;
            int r2 = 0;
            for (auto e = dual.domains.edges[v->first].cbegin(); e != dual.domains.edges[v->first].cend(); ++e) {
                for (int kc = 0; kc < dual.domains.vertices[*e].kernel.size; ++kc) {
                    if (v->second.kernel.goffset + kr <= dual.domains.vertices[*e].kernel.goffset + kc) {
                        while (nonzeros[r2] < dual.domains.vertices[*e].kernel.goffset + kc) { ++r2; }
                        GGt.vals[c] = 0;
                        for (int k = G.rows[r1]; k < G.rows[r1 + 1]; ++k) {
                            GGt.vals[c] += G.vals[k] * dFG.vals[r2 * dFG.ncols + G.cols[k]];
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

template struct TFETIConjugateSymmetric<double>;
template struct TFETIConjugateSymmetric<std::complex<double> >;

}

