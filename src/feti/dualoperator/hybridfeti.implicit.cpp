
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
: DualOperator<T>(feti), isRegularK(false)
{

}

template <typename T>
HybridFETIImplicit<T>::~HybridFETIImplicit()
{

}

template <typename T>
void HybridFETIImplicit<T>::info()
{
    if (this->infoPrinted && !feti.updated.K && !feti.updated.B) {
        return;
    }
    this->infoPrinted = true;

    DualOperatorInfo sum, min, max;
    DualOperator<T>::reduceInfo(KSolver, sum, min, max);

    eslog::info("      = ------------------------------------------------------------------------------- = \n");
    eslog::info("      = IMPLICIT HYBRID FETI OPERATOR                                                   = \n");
    DualOperator<T>::printInfo(KSolver, sum, min, max);
    eslog::info("      = ------------------------------------------------------------------------------- = \n");
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

    int nR1 = 0, nKR1 = 0;
    for (size_t di = 0; di < feti.K.size(); ++di) {
        nR1 += feti.R1[di].nrows;
        nKR1 += feti.KR1[di].nrows;
    }
    if (nR1 == 0 && nKR1 == 0) {
        eslog::error("HYBRID FETI: provide kernels for B0 gluing.\n");
    }
    isRegularK = nR1 == 0;

    _computeB0();
    _computeF0();
    _computeG0();
    _computeS0();
    g.resize(feti.cluster.gl_size);
    beta.resize(G0.nrows);
    mu.resize(feti.cluster.gl_size);

    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{B0, F0, G0, S0}\n");
        for (size_t d = 0; d < feti.B1.size(); ++d) {
            math::store(B0[d], utils::filename(utils::debugDirectory(step) + "/feti/dualop", "B0" + std::to_string(d)).c_str());
            math::store(D2C0[d], utils::filename(utils::debugDirectory(step) + "/feti/dualop", "D2C0" + std::to_string(d)).c_str());
        }
        math::store(F0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "F0").c_str());
        math::store(G0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "G0").c_str());
        math::store(S0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "S0").c_str());
    }

    _applyK(feti.f, KplusBtx);
    applyB(feti, KplusBtx, d);
    d.synchronize();
    math::add(d, T{-1}, feti.c);
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{d}\n");
        math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
    }

    info();
}

template <typename T>
void HybridFETIImplicit<T>::_apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
    }

    _applyK(Btx, KplusBtx);
    applyB(feti, KplusBtx, y);
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
        applyBt(feti, di, x, Btx[di]);
        math::copy(KplusBtx[di], feti.f[di]);
        math::add(KplusBtx[di], T{-1}, Btx[di]);
    }
   _applyK(KplusBtx, y);
}

template <typename T>
void HybridFETIImplicit<T>::BtL(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, y[di]);
    }

    _compute_beta_mu(y);

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        // y = (b - B0t * mu)
        math::spblas::applyT(y[di], T{-1}, B0[di], D2C0[di].data(), mu);
    }
}

// https://dl.acm.org/doi/pdf/10.1145/2929908.2929909

template <typename T>
void HybridFETIImplicit<T>::_applyK(std::vector<Vector_Dense<T> > &b, std::vector<Vector_Dense<T> > &x)
{
    _compute_beta_mu(b);

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        // Btx = (b - B0t * mu)
        math::spblas::applyT(b[di], T{-1}, B0[di], D2C0[di].data(), mu);
        // x = (K+)^-1 * (b - B0t * mu)
        KSolver[di].solve(b[di], x[di]);
    }

    if (!isRegularK) {
        // x += R * beta
        #pragma omp parallel for
        for (size_t di = 0; di < feti.K.size(); ++di) {
            Vector_Dense<T> _beta;
            _beta.size = origR1[di].nrows;
            _beta.vals = beta.vals + G0offset[di];
            math::blas::multiply(T{1}, origR1[di], _beta, T{1}, x[di], true);
        }
    }
}

template <typename T>
void HybridFETIImplicit<T>::_compute_beta_mu(std::vector<Vector_Dense<T> > &b)
{
    // g = B0 * (K+)^-1 * b
    math::set(g, T{0});
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Vector_Dense<T> KB0b; KB0b.resize(dKB0[di].nrows);
        math::blas::apply(KB0b, T{1}, dKB0[di], T{0}, b[di]);
        for (size_t i = 0; i < D2C0[di].size(); ++i) {
            g.vals[D2C0[di][i]] += KB0b.vals[i];
        }

        if (!isRegularK) {
            // e = Rt * b
            Vector_Dense<T> _e;
            _e.size = origR1[di].nrows;
            _e.vals = beta.vals + G0offset[di];
            math::blas::multiply(T{1}, origR1[di], b[di], T{0}, _e); // assemble -e
        }
    }

    if (!isRegularK) {
        auto &F0g = mu; // mu is tmp variable that can be used here
        F0Solver.solve(g, F0g);
        // e = G0 * F^-1 * g - e
        math::spblas::apply(beta, T{1}, G0, F0g);
        // beta = (S+)^-1 * (G0 * F0^-1 * g - e)
        Splus.solve(beta);
        // g = g - G0t * beta
        math::spblas::applyT(g, T{-1}, G0, beta);
    }
    // mu = F0^-1 * (g - G0t * beta)
    F0Solver.solve(g, mu);
}

template <typename T>
void HybridFETIImplicit<T>::_computeB0()
{
    std::vector<Matrix_Dense<T> > &R = isRegularK ? feti.KR1 : feti.R1;

    B0.clear();
    B0.resize(feti.K.size());
    D2C0.clear();
    D2C0.resize(feti.K.size());

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
    std::vector<std::vector<__ijv__> > iB0(B0.size());

    auto dual = info::mesh->domains->localDual->begin();
    std::vector<int> rows(info::mesh->clusters->size);
    for (int d1 = 0; d1 < info::mesh->domains->size; ++d1, ++dual) {
        int cluster = info::mesh->domains->cluster[d1];
        for (auto dit = dual->begin(); dit != dual->end(); ++dit) {
            if (d1 < *dit) {
                rindex.push_back(rows[cluster]);
                if (R[d1].nrows < R[*dit].nrows) {
                    rows[cluster] += R[*dit].nrows;
                } else {
                    rows[cluster] += R[d1].nrows ? R[d1].nrows : 1;
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
                        if (R[d1].nrows) {
                            for (int r = 0; r < R[d1].nrows; ++r) {
                                iB0[d1].push_back({rindex[it - dual->begin()] + r, d1index,  R[d1].vals[R[d1].ncols * r + d1index]});
                                iB0[d2].push_back({rindex[it - dual->begin()] + r, d2index, -R[d1].vals[R[d1].ncols * r + d1index]});
                            }
                        } else {
                            iB0[d1].push_back({rindex[it - dual->begin()], d1index,  (double)(d1 + 1) / info::mesh->domains->size});
                            iB0[d1].push_back({rindex[it - dual->begin()], d2index, -(double)(d1 + 1) / info::mesh->domains->size});
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int d = 0; d < info::mesh->domains->size; ++d) {
        std::sort(iB0[d].begin(), iB0[d].end());
        if (iB0[d].size()) {
            D2C0[d].push_back(iB0[d][0].i);
        }
        for (size_t i = 1; i < iB0[d].size(); ++i) {
            if (D2C0[d].back() != iB0[d][i].i) {
                D2C0[d].push_back(iB0[d][i].i);
            }
        }
        B0[d].resize(D2C0[d].size(), feti.K[d].ncols, iB0[d].size());
        B0[d].rows[0] = 0;
        B0[d].cols[0] = iB0[d][0].j;
        B0[d].vals[0] = iB0[d][0].v;
        for (size_t i = 1, r = 1; i < iB0[d].size(); ++i) {
            B0[d].cols[i] = iB0[d][i].j;
            B0[d].vals[i] = iB0[d][i].v;
            if (iB0[d][i - 1].i != iB0[d][i].i) {
                B0[d].rows[r++] = i;
            }
        }
        B0[d].rows[D2C0[d].size()] = iB0[d].size();
    }
}

template <typename T>
void HybridFETIImplicit<T>::_computeF0()
{
    if (F0.type == Matrix_Type::UNSET_INVALID_NONE) {
        std::vector<std::vector<int> > csr(feti.cluster.gl_size);
        for (size_t di = 0; di < feti.K.size(); ++di) {
            for (size_t i = 0; i < D2C0[di].size(); ++i) {
                for (size_t j = i; j < D2C0[di].size(); ++j) {
                    csr[D2C0[di][i]].push_back(D2C0[di][j]);
                }
            }
            dKB0[di].resize(B0[di].nrows, B0[di].ncols);
        }

        #pragma omp parallel for
        for (size_t i = 0; i < csr.size(); ++i) {
            utils::sortAndRemoveDuplicates(csr[i]);
        }
        int F0nnz = 0;
        for (size_t i = 0; i < csr.size(); ++i) {
            F0nnz += csr[i].size();
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

        permutationF0.resize(feti.K.size());
        for (size_t di = 0, ri = feti.cluster.gl_size; di < feti.K.size(); ++di, ++ri) {
            for (size_t i = 0; i < D2C0[di].size(); ++i) {
                for (size_t j = i; j < D2C0[di].size(); ++j) {
                    int c = std::lower_bound(csr[D2C0[di][i]].begin(), csr[D2C0[di][i]].end(), D2C0[di][j]) - csr[D2C0[di][i]].begin();
                    permutationF0[di].push_back(F0.rows[D2C0[di][i]] + c - Indexing::CSR);
                }
            }
        }

        eslog::checkpointln("FETI: SET F0");
        F0Solver.commit(F0);
        F0Solver.symbolicFactorization();
        eslog::checkpointln("FETI: F0 SYMBOLIC FACTORIZATION");
    }

    math::set(F0, T{0});
    std::vector<int> pi(feti.K.size());
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Matrix_Dense<T> dB0; dB0.resize(B0[di].nrows, B0[di].ncols);
        Matrix_Dense<T> dF0; dF0.resize(B0[di].nrows, B0[di].nrows);
        math::set(dB0, T{0});
        for (size_t i = 0; i < D2C0[di].size(); ++i) {
            for (int c = B0[di].rows[i]; c < B0[di].rows[i + 1]; ++c) {
                dB0.vals[i * dB0.ncols + B0[di].cols[c]] = B0[di].vals[c];
            }
        }
        // dF0 = B0 * (K+)^-1 * B0t // dense due to B0 is local
        KSolver[di].solve(dB0, dKB0[di]);
        math::blas::multiply(T{1}, dB0, dKB0[di], T{0}, dF0, false, true);
        // dF0 to F0
        for (size_t i = 0, k = 0; i < D2C0[di].size(); k += ++i) {
            for (size_t j = i; j < D2C0[di].size(); ++j, ++k) {
                #pragma omp atomic
                F0.vals[permutationF0[di][pi[di]++]] += dF0.vals[k];
            }
        }
    }
    eslog::checkpointln("FETI: UPDATE F0");
    F0Solver.numericalFactorization();
    eslog::checkpointln("FETI: F0 NUMERIC FACTORIZATION");
}

template <typename T>
void HybridFETIImplicit<T>::_computeG0()
{
    if (isRegularK) {
        G0.resize(0, feti.cluster.gl_size, 0);
        G0.rows[0] = 0;
        return;
    }

    int G0nnz = 0, G0rows = 0;
    for (size_t di = 0; di < feti.K.size(); ++di) {
        G0nnz += feti.R1[di].nrows * D2C0[di].size();
        G0rows += feti.R1[di].nrows;
    }
    G0.resize(G0rows, feti.cluster.gl_size, G0nnz);
    G0.rows[0] = 0;
    for (size_t di = 0, ri = 0, ci = 0; di < feti.K.size(); ++di) {
        G0offset.push_back(ri);
        for (int r = 0; r < feti.R1[di].nrows; ++r, ++ri) {
            for (size_t c = 0; c < D2C0[di].size(); ++c) {
                G0.cols[ci++] = D2C0[di][c];
            }
            G0.rows[ri + 1] = ci;
        }
    }

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        // G0 = R1 * B0t
        // dense version: math::blas::multiply(T{1}, feti.R1[di], dB0[t], T{0}, dKB0[t], false, true);
        // sparse version
        for (int r = 0, ri = G0offset[di]; r < feti.R1[di].nrows; ++r, ++ri) {
            for (size_t c = 0; c < D2C0[di].size(); ++c) {
                G0.vals[G0.rows[ri] + c] = 0;
                for (int k = B0[di].rows[c]; k < B0[di].rows[c + 1]; ++k) {
                    G0.vals[G0.rows[ri] + c] -= feti.R1[di].vals[r * feti.R1[di].ncols + B0[di].cols[k]] * B0[di].vals[k];
                }
            }
        }
    }
    eslog::checkpointln("FETI: UPDATE G0");
}

template <typename T>
void HybridFETIImplicit<T>::_computeS0()
{
    if (isRegularK) {
        S0.resize(0, 0);
        Splus.commit(S0);
        Splus.factorization();
        return;
    }

    Matrix_Dense<T> dG0, dF0G0;
    math::copy(dG0, G0);
    F0Solver.solve(dG0, dF0G0);
    S0.resize(G0.nrows, G0.nrows); S0.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    math::blas::multiply(T{1}, dG0, dF0G0, T{0}, S0, false, true); // TODO: use sparse G0?

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
    case FETIConfiguration::REGULARIZATION::ALGEBRAIC:
    case FETIConfiguration::REGULARIZATION::SVD: {
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

    eslog::checkpointln("FETI: UPDATE S0");

    Splus.commit(S0);
    Splus.factorization();
    eslog::checkpointln("FETI: S0 FACTORIZATION");
}

template class HybridFETIImplicit<double>;

}



