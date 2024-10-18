
#include "totalfeti.implicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

template <typename T>
TotalFETIImplicit<T>::TotalFETIImplicit(FETI<T> &feti)
: DualOperator<T>(feti), sparsity(DirectSparseSolver<T>::VectorSparsity::DENSE)
{

}

template <typename T>
TotalFETIImplicit<T>::~TotalFETIImplicit()
{

}

template <typename T>
void TotalFETIImplicit<T>::info()
{
    if (this->infoPrinted && !feti.updated.K && !feti.updated.B) {
        return;
    }
    this->infoPrinted = true;

    DualOperatorInfo sum, min, max;
    DualOperator<T>::reduceInfo(KSolver, sum, min, max);

    eslog::info("      = ------------------------------------------------------------------------------- = \n");
    eslog::info("      = IMPLICIT TOTAL FETI OPERATOR                                                    = \n");
    DualOperator<T>::printInfo(KSolver, sum, min, max);
    if (feti.configuration.exhaustive_info > 2) {
        eslog::info("      =  A = K, B = REG(K), Y = A+                                                      = \n");

        std::vector<Matrix_Dense<T> > A(feti.K.size()), B(feti.K.size()), Y(feti.K.size());
        std::vector<T> min(6, 1e99), max(6), sum(6), normK(feti.K.size());

        for (size_t di = 0; di < feti.K.size(); ++di) {
            math::copy(A[di], Kplus[di]);
            math::copy(B[di], feti.K[di]);

            Matrix_Dense<T> I; I.resize(feti.K[di].nrows, feti.K[di].ncols);
            math::eye(I, T{1});
            KSolver[di].solve(I, Y[di]);
        }

        { // || DIAG(A) ||
            for (size_t di = 0; di < feti.K.size(); ++di) {
                for (int r = 0; r < feti.K[di].nrows; ++r) {
                    normK[di] += feti.K[di].vals[feti.K[di].rows[r - Indexing::CSR]] * feti.K[di].vals[feti.K[di].rows[r - Indexing::CSR]];
                }
                normK[di] = std::sqrt(normK[di]);
                min[0] = std::min(min[0], normK[di]);
                max[0] = std::max(max[0], normK[di]);
                sum[0] += normK[di];
            }
        }

        { // || BY - I ||
            Matrix_Dense<T> BY, I;
            for (size_t di = 0; di < feti.K.size(); ++di) {
                BY.resize(B[di].nrows, B[di].ncols);
                I.resize(B[di].nrows, B[di].ncols);
                math::eye(I, T{1});
                math::blas::multiply(T{1}, A[di], Y[di], T{0}, BY);
                math::add(BY, T{-1}, I);

                T norm = math::norm(BY);
                min[1] = std::min(min[1], norm);
                max[1] = std::max(max[1], norm);
                sum[1] += norm;
            }
        }
        { // || AYA - A ||
            Matrix_Dense<T> AY, AYA;
            for (size_t di = 0; di < feti.K.size(); ++di) {
                AY.resize(A[di].nrows, A[di].ncols);
                AYA.resize(A[di].nrows, A[di].ncols);
                math::blas::multiply(T{1}, A[di], Y[di], T{0}, AY);
                math::blas::multiply(T{1}, AY, A[di], T{0}, AYA);
                math::add(AYA, T{-1}, A[di]);

                T norm = math::norm(AYA);
                min[2] = std::min(min[2], norm);
                max[2] = std::max(max[2], norm);
                sum[2] += norm;
            }
        }
        { // || YAY - Y ||
            Matrix_Dense<T> YA, YAY;
            for (size_t di = 0; di < feti.K.size(); ++di) {
                YA.resize(A[di].nrows, A[di].ncols);
                YAY.resize(A[di].nrows, A[di].ncols);
                math::blas::multiply(T{1}, Y[di], A[di], T{0}, YA);
                math::blas::multiply(T{1}, YA, Y[di], T{0}, YAY);
                math::add(YAY, T{-1}, Y[di]);

                T norm = math::norm(YAY);
                min[3] = std::min(min[3], norm);
                max[3] = std::max(max[3], norm);
                sum[3] += norm;
            }
        }
        { // || AY - (AY)' ||
            Matrix_Dense<T> AY, AYT;
            for (size_t di = 0; di < feti.K.size(); ++di) {
                AY.resize(A[di].nrows, A[di].ncols);
                AYT.resize(A[di].nrows, A[di].ncols);
                math::blas::multiply(T{1}, A[di], Y[di], T{0}, AY);
                math::blas::multiply(T{1}, A[di], Y[di], T{0}, AYT, true, true);
                math::add(AY, T{-1}, AYT);

                T norm = math::norm(AY);
                min[4] = std::min(min[4], norm);
                max[4] = std::max(max[4], norm);
                sum[4] += norm;
            }
        }
        { // || YA - (YA)' ||
            Matrix_Dense<T> YA, YAT;
            for (size_t di = 0; di < feti.K.size(); ++di) {
                YA.resize(A[di].nrows, A[di].ncols);
                YAT.resize(A[di].nrows, A[di].ncols);
                math::blas::multiply(T{1}, Y[di], A[di], T{0}, YA);
                math::blas::multiply(T{1}, Y[di], A[di], T{0}, YAT, true, true);
                math::add(YA, T{-1}, YAT);

                T norm = math::norm(YA);
                min[5] = std::min(min[5], norm);
                max[5] = std::max(max[5], norm);
                sum[5] += norm;
            }
        }
        Communication::allReduce(min.data(), nullptr, 6, MPITools::getType<T>().mpitype, MPI_MIN);
        Communication::allReduce(max.data(), nullptr, 6, MPITools::getType<T>().mpitype, MPI_MAX);
        Communication::allReduce(sum.data(), nullptr, 6, MPITools::getType<T>().mpitype, MPI_SUM);
        eslog::info("      =  || A ||                                        %.2e <%.2e - %.2e>  = \n", (double)sum[0] / feti.sinfo.domains, min[0], max[0]);
        eslog::info("      =  || BY - I ||                                   %.2e <%.2e - %.2e>  = \n", (double)sum[1] / feti.sinfo.domains, min[1], max[1]);
        eslog::info("      =  || AYA - A ||                                  %.2e <%.2e - %.2e>  = \n", (double)sum[2] / feti.sinfo.domains, min[2], max[2]);
        eslog::info("      =  || YAY - Y ||                                  %.2e <%.2e - %.2e>  = \n", (double)sum[3] / feti.sinfo.domains, min[3], max[3]);
        eslog::info("      =  || AY - (AY)' ||                               %.2e <%.2e - %.2e>  = \n", (double)sum[4] / feti.sinfo.domains, min[4], max[4]);
        eslog::info("      =  || YA - (YA)' ||                               %.2e <%.2e - %.2e>  = \n", (double)sum[5] / feti.sinfo.domains, min[5], max[5]);
    }
    eslog::info("      = ------------------------------------------------------------------------------- = \n");
}

/*
 * prepare buffers and call symbolic factorization that is independent on the Kplus values
 */
template <typename T>
void TotalFETIImplicit<T>::set(const step::Step &step)
{
    sparsity = feti.configuration.partial_dual ? DirectSparseSolver<T>::VectorSparsity::SPARSE_RHS | DirectSparseSolver<T>::VectorSparsity::SPARSE_SOLUTION : DirectSparseSolver<T>::VectorSparsity::DENSE;

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

        int suffix = 0;
        if (sparsity != DirectSparseSolver<T>::VectorSparsity::DENSE) {
            suffix = *std::min_element(feti.B1[di].cols, feti.B1[di].cols + feti.B1[di].nnz);
        }

        KSolver[di].symbolicFactorization(suffix);
    }
    eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");
}

template <typename T>
void TotalFETIImplicit<T>::update(const step::Step &step)
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

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        KSolver[di].solve(feti.f[di], KplusBtx[di], sparsity);
    }
    applyB(feti, KplusBtx, d);
    math::add(d, T{-1}, feti.c);
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{d}\n");
        math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
    }

    info();
}


template <typename T>
void TotalFETIImplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
        KSolver[di].solve(Btx[di], KplusBtx[di], sparsity);
    }
    applyB(feti, KplusBtx, y);
}

template <typename T>
void TotalFETIImplicit<T>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
        math::copy(KplusBtx[di], feti.f[di]);
        math::add(KplusBtx[di], T{-1}, Btx[di]);
        KSolver[di].solve(KplusBtx[di], y[di]);
    }
}

template class TotalFETIImplicit<double>;

}
