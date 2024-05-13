
#include "hybridfeti.implicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
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
}

template <typename T>
void HybridFETIImplicit<T>::update(const step::Step &step)
{
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
    Matrix_Dense<T> S0; S0.resize(G0.nrows, G0.nrows);
    math::blas::multiply(T{1}, dG0, dF0G0, T{0}, S0, false, true); // TODO: use sparse G0?

    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{S0}\n");
        math::store(S0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "S0").c_str());
    }

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

    // make S0 regular
    for (size_t di = 0, ri = 0; di < feti.K.size(); ++di) {
        for (int r = 0; r < feti.R1[di].nrows; ++r, ++ri) {
            for (size_t dj = 0; dj < feti.K.size(); ++dj) {
                if (r < feti.R1[dj].nrows) {
                    S0.vals[ri * S0.ncols + G0offset[dj] + r] += 0.5 * avg;
                }
            }
        }
    }

    eslog::checkpointln("FETI: S0 ASSEMBLE");

    Splus.commit(std::move(S0));
    Splus.factorization();

    eslog::checkpointln("FETI: S0 FACTORIZATION");

    _applyK(feti.f, KplusBtx);
    applyB(feti, KplusBtx, d);
    math::add(d, T{-1}, feti.c);
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{d}\n");
        math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
    }
}

template <typename T>
void HybridFETIImplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
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
        _e.size = feti.R1[di].nrows;
        _e.vals = beta.vals + G0offset[di];
        math::blas::multiply(T{1}, feti.R1[di], x[di], T{0}, _e); // assemble -e
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
        _beta.size = feti.R1[di].nrows;
        _beta.vals = beta.vals + G0offset[di];
        math::blas::multiply(T{1}, feti.R1[di], _beta, T{1}, y[di], true);
    }
}

template class HybridFETIImplicit<double>;

}



