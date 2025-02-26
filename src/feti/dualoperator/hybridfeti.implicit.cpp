
#include "hybridfeti.implicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

template <typename T>
HybridFETIImplicit<T>::HybridFETIImplicit(FETI<T> &feti)
: DualOperator<T>(feti)
{

}

template <typename T>
HybridFETIImplicit<T>::~HybridFETIImplicit()
{

}

template <typename T>
void HybridFETIImplicit<T>::info()
{
    DualOperatorInfo sum, min, max;
    DualOperator<T>::reduceInfo(KSolver, sum, min, max);

    eslog::info(" = IMPLICIT HYBRID FETI OPERATOR                                                             = \n");
    DualOperator<T>::printInfo(KSolver, sum, min, max);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

/*
 * prepare buffers and call symbolic factorization that is independent on the Kplus values
 */
template <typename T>
void HybridFETIImplicit<T>::set(const step::Step &step)
{
    Kplus.resize(feti.K.size());
    Btx.resize(feti.K.size());
    KplusBtx.resize(feti.K.size());
    KSolver.resize(feti.K.size());
    dKB0.resize(feti.K.size());

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Kplus[di].type = feti.K[di].type;
        Kplus[di].shape = feti.K[di].shape;
        math::combine(Kplus[di], feti.K[di], feti.RegMat[di]);
        Btx[di].resize(feti.K[di].nrows);
        KplusBtx[di].resize(feti.K[di].nrows);
        math::set(Btx[di], T{0});
    }

    eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        KSolver[di].commit(Kplus[di]);
        KSolver[di].symbolicFactorization();
    }
    eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");
}

template <typename T>
void HybridFETIImplicit<T>::update(const step::Step &step)
{
    if (feti.updated.B) {
        d.resize();
    }

    if (feti.updated.K) {
        _computeB0();
    }

    int dB0max = 0, dKmax = 0;
    std::vector<std::vector<int> > csr(feti.cluster.gl_size);
    for (size_t di = 0; di < feti.K.size(); ++di) {
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            for (size_t j = i; j < feti.D2C0[di].size(); ++j) {
                csr[feti.D2C0[di][i]].push_back(feti.D2C0[di][j]);
            }
        }
        dB0max = std::max(dB0max, feti.B0[di].nrows);
        dKmax = std::max(dKmax, feti.K[di].nrows);
        dKB0[di].resize(feti.B0[di].nrows, feti.B0[di].ncols);
    }

    dB0.resize(info::env::threads);
    dF0.resize(info::env::threads);
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        dB0[t].resize(dB0max, dKmax);
        dF0[t].resize(dB0max, dB0max);
        math::set(dB0[t], T{0});
        math::set(dF0[t], T{0});
    }

    #pragma omp parallel for
    for (size_t i = 0; i < csr.size(); ++i) {
        utils::sortAndRemoveDuplicates(csr[i]);
    }
    int F0nnz = 0, G0nnz = 0, G0rows = 0;
    for (size_t i = 0; i < csr.size(); ++i) {
        F0nnz += csr[i].size();
    }
    for (size_t di = 0; di < feti.K.size(); ++di) {
        G0nnz += feti.R1[di].nrows * feti.D2C0[di].size();
        G0rows += feti.R1[di].nrows;
    }
    F0.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    F0.shape = Matrix_Shape::UPPER;
    F0.resize(feti.cluster.gl_size, feti.cluster.gl_size, F0nnz);
    F0.rows[0] = Indexing::CSR;
    for (size_t i = 0; i < csr.size(); ++i) {
        F0.rows[i + 1] = F0.rows[i] + std::max(csr[i].size(), 1UL);
        for (size_t j = 0; j < csr[i].size(); ++j) {
            F0.cols[F0.rows[i] - Indexing::CSR + j] = csr[i][j] + Indexing::CSR;
        }
    }
    G0.resize(G0rows, feti.cluster.gl_size, G0nnz);
    G0.rows[0] = 0;
    for (size_t di = 0, ri = 0, ci = 0; di < feti.K.size(); ++di) {
        G0offset.push_back(ri);
        for (int r = 0; r < feti.R1[di].nrows; ++r, ++ri) {
            for (size_t c = 0; c < feti.D2C0[di].size(); ++c) {
                G0.cols[ci++] = feti.D2C0[di][c];
            }
            G0.rows[ri + 1] = ci;
        }
    }

    permutation.resize(feti.K.size());
    for (size_t di = 0, ri = feti.cluster.gl_size; di < feti.K.size(); ++di, ++ri) {
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            for (size_t j = i; j < feti.D2C0[di].size(); ++j) {
                int c = std::lower_bound(csr[feti.D2C0[di][i]].begin(), csr[feti.D2C0[di][i]].end(), feti.D2C0[di][j]) - csr[feti.D2C0[di][i]].begin();
                permutation[di].push_back(F0.rows[feti.D2C0[di][i]] + c - Indexing::CSR);
            }
        }
    }

    g.resize(G0.ncols);
    beta.resize(G0.nrows);
    mu.resize(G0.ncols);
    eslog::checkpointln("FETI: SET F0, G0");

    F0Solver.commit(F0);
    F0Solver.symbolicFactorization();
    eslog::checkpointln("FETI: F0 SYMBOLIC FACTORIZATION");

    if (feti.updated.K) {
        #pragma omp parallel for
        for (size_t di = 0; di < feti.K.size(); ++di) {
            math::sumCombined(Kplus[di], T{1}, feti.K[di], feti.RegMat[di]);
        }
        eslog::checkpointln("FETI: UPDATE TOTAL-FETI OPERATOR");
    }
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{Kplus}\n");
        for (size_t di = 0; di < feti.K.size(); ++di) {
            math::store(Kplus[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(di)).c_str()).c_str());
        }
    }

    if (feti.updated.K) {
        #pragma omp parallel for
        for (size_t di = 0; di < feti.K.size(); ++di) {
            KSolver[di].numericalFactorization();
        }
        eslog::checkpointln("FETI: TFETI NUMERICAL FACTORIZATION");
    }

    math::set(F0, T{0});
    std::vector<int> pi(feti.K.size());
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        int t = omp_get_thread_num();
        dB0[t].resize(feti.B0[di].nrows, feti.B0[di].ncols);
        dF0[t].resize(feti.B0[di].nrows, feti.B0[di].nrows);
        math::set(dB0[t], T{0});
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            for (int c = feti.B0[di].rows[i]; c < feti.B0[di].rows[i + 1]; ++c) {
                dB0[t].vals[i * dB0[t].ncols + feti.B0[di].cols[c]] = feti.B0[di].vals[c];
            }
        }
        // dF0 = B0 * (K+)^-1 * B0t // dense due to B0 is local
        KSolver[di].solve(dB0[t], dKB0[di]);
        math::blas::multiply(T{1}, dB0[t], dKB0[di], T{0}, dF0[t], false, true);
        // dF0 to F0
        for (size_t i = 0, k = 0; i < feti.D2C0[di].size(); k += ++i) {
            for (size_t j = i; j < feti.D2C0[di].size(); ++j, ++k) {
                #pragma omp atomic
                F0.vals[permutation[di][pi[di]++]] += dF0[t].vals[k];
            }
        }

        // G0 = R1 * B0t
        // dense version: math::blas::multiply(T{1}, feti.R1[di], dB0[t], T{0}, dKB0[t], false, true);
        // sparse version
        for (int r = 0, ri = G0offset[di]; r < feti.R1[di].nrows; ++r, ++ri) {
            for (size_t c = 0; c < feti.D2C0[di].size(); ++c) {
                G0.vals[G0.rows[ri] + c] = 0;
                for (int k = feti.B0[di].rows[c]; k < feti.B0[di].rows[c + 1]; ++k) {
                    G0.vals[G0.rows[ri] + c] -= feti.R1[di].vals[r * feti.R1[di].ncols + feti.B0[di].cols[k]] * feti.B0[di].vals[k];
                }
            }
        }
    }
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{F0, G0}\n");
        math::store(F0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "F0").c_str());
        math::store(G0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "G0").c_str());
    }
    eslog::checkpointln("FETI: UPDATE F0, G0");

    F0Solver.numericalFactorization();
    eslog::checkpointln("FETI: F0 NUMERICAL FACTORIZATION");

    Matrix_Dense<T> dG0, dF0G0;
    math::copy(dG0, G0);
    F0Solver.solve(dG0, dF0G0);
    Matrix_Dense<T> S0; S0.resize(G0.nrows, G0.nrows); S0.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    math::blas::multiply(T{1}, dG0, dF0G0, T{0}, S0, false, true); // TODO: use sparse G0?

    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{S0}\n");
        math::store(S0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "S0").c_str());
    }

    origR1.resize(feti.K.size());
    // make S0 regular
    switch (feti.configuration.regularization) {
    case FETIConfiguration::REGULARIZATION::ANALYTIC: {
        T avg = 0;
        for (int r = 0; r < S0.nrows; ++r) {
            avg += S0.vals[S0.nrows * r + r];
        }
        avg /= S0.nrows;

        // kernels: [3, 2, 3]
        // 1 0 0   *  1 0 0 1 0 1 0 0  =   1 0 0 1 0 1 0 0
        // 0 1 0      0 1 0 0 1 0 1 0      0 1 0 0 1 0 1 0
        // 0 0 1      0 0 1 0 0 0 0 1      0 0 1 0 0 0 0 1
        // 1 0 0                           1 0 0 1 0 1 0 0
        // 0 1 0                           0 1 0 0 1 0 1 0
        // 1 0 0                           1 0 0 1 0 1 0 0
        // 0 1 0                           0 1 0 0 1 0 1 0
        // 0 0 1                           0 0 1 0 0 0 0 1

        for (size_t di = 0, ri = 0; di < feti.K.size(); ++di) {
            origR1[di].shallowCopy(feti.R1[di]);
            for (int r = 0; r < feti.R1[di].nrows; ++r, ++ri) {
                for (size_t dj = 0; dj < feti.K.size(); ++dj) {
                    if (r < feti.R1[dj].nrows) {
                        S0.vals[ri * S0.ncols + G0offset[dj] + r] += 0.5 * avg;
                    }
                }
            }
        }
    } break;
    case FETIConfiguration::REGULARIZATION::ALGEBRAIC: {
        int maxDefect = 0;
        for (size_t di = 0; di < feti.K.size(); ++di) {
            maxDefect = std::max(maxDefect, feti.R1[di].nrows);
            if (maxDefect != feti.R1[di].nrows) {
                eslog::error("Some domains are regular\n");
            }
        }
        Matrix_Dense<T> R;
        Matrix_IJV<T> regMat;
        math::getKernel<T, int>(S0, R, regMat, maxDefect, feti.configuration.sc_size);
        for (int i = 0; i < regMat.nnz; ++i) {
            S0.vals[(regMat.rows[i] - Indexing::IJV) * S0.ncols + regMat.cols[i] - Indexing::IJV] += regMat.vals[i];
        }
        for (size_t di = 0; di < feti.K.size(); ++di) {
            if (feti.configuration.projector_opt & FETIConfiguration::PROJECTOR_OPT::FULL) {
                origR1[di].shallowCopy(feti.R1[di]);
            } else {
                // orthogonalize R for HFETI projector
                origR1[di].resize(feti.R1[di].nrows, feti.R1[di].ncols);
                math::copy(origR1[di], feti.R1[di]);
                Matrix_Dense<T> _R;
                math::lapack::submatrix<T, int>(R, _R, 0, R.nrows, di * R.nrows, (di + 1) * R.nrows);
                math::blas::multiply(T{1}, _R, origR1[di], T{0}, feti.R1[di]);
            }
        }
    } break;
    }

    eslog::checkpointln("FETI: S0 ASSEMBLE");

    Splus.commit(std::move(S0));
    Splus.factorization();

    eslog::checkpointln("FETI: S0 FACTORIZATION");

    _applyK(feti.f, KplusBtx);
    applyB(feti, KplusBtx, d);
    d.synchronize();
    math::add(d, T{-1}, feti.c);
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{d}\n");
        math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
    }
}

template <typename T>
void HybridFETIImplicit<T>::_computeB0()
{
    feti.B0.clear();
    feti.D2C0.clear();
    feti.B0.resize(feti.K.size());
    feti.D2C0.resize(feti.K.size());

    struct __ijv__ {
        int i,j; double v;

        bool operator<(const __ijv__ &other) const {
            if (i == other.i) {
                return j < other.j;
            }
            return i < other.i;
        }
    };

    std::vector<int> rindex;
    std::vector<std::vector<__ijv__> > B0(feti.B0.size());

    auto dual = info::mesh->domains->localDual->begin();
    std::vector<int> rows(info::mesh->clusters->size);
    for (int d1 = 0; d1 < info::mesh->domains->size; ++d1, ++dual) {
        int cluster = info::mesh->domains->cluster[d1];
        for (auto dit = dual->begin(); dit != dual->end(); ++dit) {
            if (d1 < *dit) {
                rindex.push_back(rows[cluster]);
                if (feti.R1[d1].nrows < feti.R1[*dit].nrows) {
                    rows[cluster] += feti.R1[*dit].nrows;
                } else {
                    rows[cluster] += feti.R1[d1].nrows ? feti.R1[d1].nrows : 1;
                }
            } else {
                auto dualbegin = info::mesh->domains->localDual->begin();
                auto it = std::lower_bound((dualbegin + *dit)->begin(), (dualbegin + *dit)->end(), d1);
                rindex.push_back(rindex[it - dualbegin->begin()]);
            }
        }
    }
    feti.cluster.gl_size = rows.front(); // TODO: more clusters

    dual = info::mesh->domains->localDual->begin();
    for (auto dmap = feti.decomposition->dmap->cbegin(); dmap != feti.decomposition->dmap->cend(); ++dmap) {
        for (auto di1 = dmap->begin(); di1 != dmap->end(); ++di1) {
            for (auto di2 = di1 + 1; di2 != dmap->end(); ++di2) {
                if (feti.decomposition->ismy(di1->domain) && feti.decomposition->ismy(di2->domain)) {
                    auto it = std::lower_bound(
                            (dual + (di1->domain - feti.decomposition->dbegin))->begin(),
                            (dual + (di1->domain - feti.decomposition->dbegin))->end(),
                            di2->domain - feti.decomposition->dbegin);

                    if (it != (dual + (di1->domain - feti.decomposition->dbegin))->end() && *it == di2->domain - feti.decomposition->dbegin) {
                        int d1, d2, d1index, d2index;
                        if (di1->domain < di2->domain) {
                            d1 = di1->domain - feti.decomposition->dbegin;
                            d2 = di2->domain - feti.decomposition->dbegin;
                            d1index = di1->index;
                            d2index = di2->index;
                        } else {
                            d1 = di2->domain - feti.decomposition->dbegin;
                            d2 = di1->domain - feti.decomposition->dbegin;
                            d1index = di2->index;
                            d2index = di1->index;
                        }
                        if (feti.R1[d1].nrows) {
                            for (int r = 0; r < feti.R1[d1].nrows; ++r) {
                                B0[d1].push_back({rindex[it - dual->begin()] + r, d1index,  feti.R1[d1].vals[feti.R1[d1].ncols * r + d1index]});
                                B0[d2].push_back({rindex[it - dual->begin()] + r, d2index, -feti.R1[d1].vals[feti.R1[d1].ncols * r + d1index]});
                            }
                        } else {
                            B0[d1].push_back({rindex[it - dual->begin()], d1index,  (double)(d1 + 1) / info::mesh->domains->size});
                            B0[d1].push_back({rindex[it - dual->begin()], d2index, -(double)(d1 + 1) / info::mesh->domains->size});
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int d = 0; d < info::mesh->domains->size; ++d) {
        std::sort(B0[d].begin(), B0[d].end());
        if (B0[d].size()) {
            feti.D2C0[d].push_back(B0[d][0].i);
        }
        for (size_t i = 1; i < B0[d].size(); ++i) {
            if (feti.D2C0[d].back() != B0[d][i].i) {
                feti.D2C0[d].push_back(B0[d][i].i);
            }
        }
        feti.B0[d].resize(feti.D2C0[d].size(), feti.K[d].ncols, B0[d].size());
        feti.B0[d].rows[0] = 0;
        feti.B0[d].cols[0] = B0[d][0].j;
        feti.B0[d].vals[0] = B0[d][0].v;
        for (size_t i = 1, r = 1; i < B0[d].size(); ++i) {
            feti.B0[d].cols[i] = B0[d][i].j;
            feti.B0[d].vals[i] = B0[d][i].v;
            if (B0[d][i - 1].i != B0[d][i].i) {
                feti.B0[d].rows[r++] = i;
            }
        }
        feti.B0[d].rows[feti.D2C0[d].size()] = B0[d].size();
    }
}

template <typename T>
void HybridFETIImplicit<T>::_apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    // FETI
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        math::set(Btx[di], T{0});
        math::spblas::applyT(Btx[di], T{1}, feti.B1[di], feti.D2C[di].data(), x);
    }

    _applyK(Btx, KplusBtx);

    // y += FETI + HFETI
    math::set(y, T{0});
    for (size_t di = 0; di < feti.K.size(); ++di) {
        math::spblas::apply(y, T{1}, feti.B1[di], feti.D2C[di].data(), KplusBtx[di]);
    }
}

template <typename T>
void HybridFETIImplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}

template <typename T>
void HybridFETIImplicit<T>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;
    for (int r = 0; r < x.nrows; ++r) {
        _x.vals = x.vals + x.ncols * r;
        _y.vals = y.vals + y.ncols * r;
        _apply(_x, _y);
    }
    y.synchronize();
}

template <typename T>
void HybridFETIImplicit<T>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        math::copy(Btx[di], feti.f[di]);
        math::spblas::applyT(Btx[di], T{-1}, feti.B1[di], feti.D2C[di].data(), x);
    }

   _applyK(Btx, y);
}

template <typename T>
void HybridFETIImplicit<T>::_applyK(std::vector<Vector_Dense<T> > &x, std::vector<Vector_Dense<T> > &y)
{
    // HYBRID FETI
    math::set(g, T{0});
    for (size_t di = 0; di < feti.K.size(); ++di) {
        // g = B0 * (K+)^-1 * B1t * x
        math::blas::apply(y[di], T{1}, dKB0[di], T{0}, x[di]);
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            g.vals[feti.D2C0[di][i]] += y[di].vals[i];
        }
        // e = Rt * b
        Vector_Dense<T> _e;
        _e.size = origR1[di].nrows;
        _e.vals = beta.vals + G0offset[di];
        math::blas::multiply(T{1}, origR1[di], x[di], T{0}, _e); // assemble -e
    }

    auto &F0g = mu; // mu is tmp variable that can be used here
    F0Solver.solve(g, F0g);
    // e = G0 * F^-1 * g - e
    math::spblas::apply(beta, T{1}, G0, F0g);
    // beta = (S+)^-1 * (G0 * F0^-1 * g - e)
    Splus.solve(beta);
    // g = g - G0t * beta
    math::spblas::applyT(g, T{-1}, G0, beta);
    // mu = F0^-1 * (g - G0t * beta)
    F0Solver.solve(g, mu);

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        // Btx = (B1t * x - B0t * mu)
        math::spblas::applyT(x[di], T{-1}, feti.B0[di], feti.D2C0[di].data(), mu);
        // x = (K+)^-1 * (B1t * x - B0t * mu)
        KSolver[di].solve(x[di], y[di]);
        // x += R * beta
        Vector_Dense<T> _beta;
        _beta.size = origR1[di].nrows;
        _beta.vals = beta.vals + G0offset[di];
        math::blas::multiply(T{1}, origR1[di], _beta, T{1}, y[di], true);
    }
}

template class HybridFETIImplicit<double>;

}



