
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
    std::vector<std::vector<int> > csr(feti.cluster.gl_size + feti.sinfo.R1size);
    for (size_t di = 0, ri = feti.cluster.gl_size; di < feti.K.size(); ++di, ++ri) {
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            for (size_t j = i; j < feti.D2C0[di].size(); ++j) {
                csr[feti.D2C0[di][i]].push_back(feti.D2C0[di][j]);
            }
            for (int j = 0; j < feti.R1[di].nrows; ++j) {
                csr[feti.D2C0[di][i]].push_back(ri + j);
            }
        }
        dB0max = std::max(dB0max, feti.B0[di].nrows);
        dKmax = std::max(dKmax, feti.K[di].nrows);
    }

    dB0.resize(info::env::threads);
    dKB0.resize(info::env::threads);
    dF0.resize(info::env::threads);
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        dB0[t].resize(dB0max, dKmax);
        dKB0[t].resize(dB0max, dKmax);
        dF0[t].resize(dB0max, dB0max);
        math::set(dB0[t], T{0});
        math::set(dF0[t], T{0});
    }

    #pragma omp parallel for
    for (size_t i = 0; i < csr.size(); ++i) {
        utils::sortAndRemoveDuplicates(csr[i]);
    }
    int nnz = 0;
    for (size_t i = 0; i < csr.size(); ++i) {
        nnz += std::max(csr[i].size(), 1UL);
    }
    HFETIDual.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    HFETIDual.resize(csr.size(), csr.size(), nnz);
    HFETIDual.rows[0] = Indexing::CSR;
    for (size_t i = 0; i < csr.size(); ++i) {
        HFETIDual.rows[i + 1] = HFETIDual.rows[i] + std::max(csr[i].size(), 1UL);
        for (size_t j = 0; j < csr[i].size(); ++j) {
            HFETIDual.cols[HFETIDual.rows[i] - Indexing::CSR + j] = csr[i][j] + Indexing::CSR;
        }
        if (csr[i].size() == 0UL) {
            HFETIDual.cols[HFETIDual.rows[i] - Indexing::CSR] = i + Indexing::CSR;
        }
    }

    permutation.resize(feti.K.size());
    for (size_t di = 0, ri = feti.cluster.gl_size; di < feti.K.size(); ++di, ++ri) {
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            for (size_t j = i; j < feti.D2C0[di].size(); ++j) {
                int c = std::lower_bound(csr[feti.D2C0[di][i]].begin(), csr[feti.D2C0[di][i]].end(), feti.D2C0[di][j]) - csr[feti.D2C0[di][i]].begin();
                permutation[di].push_back(HFETIDual.rows[feti.D2C0[di][i]] + c - Indexing::CSR);
            }
            for (int j = 0; j < feti.R1[di].nrows; ++j) {
                int c = std::lower_bound(csr[feti.D2C0[di][i]].begin(), csr[feti.D2C0[di][i]].end(), ri + j) - csr[feti.D2C0[di][i]].begin();
                permutation[di].push_back(HFETIDual.rows[feti.D2C0[di][i]] + c - Indexing::CSR);
            }
        }
    }

    HFETISolver.commit(HFETIDual);
    HFETISolver.symbolicFactorization();
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

    math::set(HFETIDual, T{0});
    std::vector<int> pi(feti.K.size());
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        int t = omp_get_thread_num();
        dB0[t].resize(feti.B0[di].nrows, feti.B0[di].ncols);
        dKB0[t].resize(feti.B0[di].nrows, feti.B0[di].ncols);
        dF0[t].resize(feti.B0[di].nrows, feti.B0[di].nrows);
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            for (int c = feti.B0[di].rows[i]; c < feti.B0[di].rows[i + 1]; ++c) {
                dB0[t].vals[i * dB0[t].ncols + feti.B0[di].cols[c]] = feti.B0[di].vals[c];
            }
        }
        KSolver[di].solve(dB0[t], dKB0[t]);
        math::blas::multiply(T{1}, dB0[t], dKB0[t], T{0}, dF0[t], false, true);
        math::store(dKB0[t], (std::string("dKB0") + std::to_string(di)).c_str());
        math::blas::multiply(T{1}, feti.R1[di], dB0[t], T{0}, dKB0[t], false, true);
        math::store(dKB0[t], (std::string("dRB0") + std::to_string(di)).c_str());
        for (size_t i = 0, k = 0; i < feti.D2C0[di].size(); k += ++i) {
            for (size_t j = i; j < feti.D2C0[di].size(); ++j, ++k) {
                HFETIDual.vals[permutation[di][pi[di]++]] += dF0[t].vals[k];
            }
            for (int j = 0; j < feti.R1[di].nrows; ++j) {
                printf("add %lu %f\n", di, -dKB0[t].vals[j * dKB0[t].ncols + i]);
                HFETIDual.vals[permutation[di][pi[di]++]] += -dKB0[t].vals[j * dKB0[t].ncols + i];
            }
        }

        math::store(dB0[t], (std::string("dB0") + std::to_string(di)).c_str());
        math::store(dF0[t], (std::string("dF0") + std::to_string(di)).c_str());
        for (size_t i = 0; i < feti.D2C0[di].size(); ++i) {
            for (int c = feti.B0[di].rows[i]; c < feti.B0[di].rows[i + 1]; ++c) {
                dB0[t].vals[i * dB0[t].ncols + feti.B0[di].cols[c]] = 0;
            }
        }

        KSolver[di].solve(feti.f[di], KplusBtx[di]);
    }
    math::store(HFETIDual, "hfeti");
    HFETISolver.numericalFactorization();

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
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
        KSolver[di].solve(Btx[di], KplusBtx[di]);
    }
    applyB(feti, KplusBtx, y);
}


template <typename T>
void HybridFETIImplicit<T>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
        math::copy(KplusBtx[di], feti.f[di]);
        math::add(KplusBtx[di], T{-1}, Btx[di]);
        KSolver[di].solve(KplusBtx[di], y[di]);
    }
}

template class HybridFETIImplicit<double>;

}



