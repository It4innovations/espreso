
#include "dualoperator.h"
#include "totalfeti.implicit.h"
#include "totalfeti.explicit.h"
#include "totalfeti.explicit.sc.h"
#include "totalfeti.gpu.h"
#include "hybridfeti.implicit.h"
#include "feti/projector/projector.h"

#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "wrappers/mpi/communication.h"

#include <climits>

namespace espreso {

template<typename T>
DualOperator<T>* DualOperator<T>::create(FETI<T> &feti, const step::Step &step)
{
    DualOperator<T>* dual = nullptr;
    switch (feti.configuration.method) {
    case FETIConfiguration::METHOD::TOTAL_FETI:
        switch (feti.configuration.dual_operator) {
        case FETIConfiguration::DUAL_OPERATOR::IMPLICIT:
            eslog::info(" = DUAL OPERATOR                                                         IMPLICIT TOTAL FETI = \n");
            dual = new TotalFETIImplicit<T>(feti);
            break;
        case FETIConfiguration::DUAL_OPERATOR::EXPLICIT:
            eslog::info(" = DUAL OPERATOR                                                         EXPLICIT TOTAL FETI = \n");
            dual = new TotalFETIExplicit<T>(feti);
            break;
        case FETIConfiguration::DUAL_OPERATOR::EXPLICIT_GPU:
            if (!gpu::mgm::is_linked()) {
                eslog::globalerror("GPU acceleration is not supported: GPU support is not built.\n");
            }
            if (!DirectSparseSolver<T>::provideFactors()) {
                eslog::globalerror("GPU acceleration is not supported: Third party sparse solver does not provide factors.\n");
            }
            eslog::info(" = DUAL OPERATOR                                                  EXPLICIT TOTAL FETI ON GPU = \n");
            dual = new TotalFETIGpu<T,int>(feti, DualOperatorStrategy::EXPLICIT);
            break;
        case FETIConfiguration::DUAL_OPERATOR::IMPLICIT_GPU:
            if (!gpu::mgm::is_linked()) {
                eslog::globalerror("GPU acceleration is not supported: GPU support is not built.\n");
            }
            if (!DirectSparseSolver<T>::provideFactors()) {
                eslog::globalerror("GPU acceleration is not supported: Third party sparse solver does not provide factors.\n");
            }
            eslog::info(" = DUAL OPERATOR                                                  IMPLICIT TOTAL FETI ON GPU = \n");
            dual = new TotalFETIGpu<T,int>(feti, DualOperatorStrategy::IMPLICIT);
            break;
        case FETIConfiguration::DUAL_OPERATOR::EXPLICIT_SC:
            eslog::info(" = DUAL OPERATOR                                  EXPLICIT TOTAL FETI USING SCHUR COMPLEMENT = \n");
            dual = new TotalFETIExplicitSc<T,int>(feti);
            break;
        }
        break;

    case FETIConfiguration::METHOD::HYBRID_FETI:
        switch (feti.configuration.dual_operator) {
        case FETIConfiguration::DUAL_OPERATOR::IMPLICIT:
            eslog::info(" = DUAL OPERATOR                                                        IMPLICIT HYBRID FETI = \n");
            dual = new HybridFETIImplicit<T>(feti);
            break;
        case FETIConfiguration::DUAL_OPERATOR::EXPLICIT:
        case FETIConfiguration::DUAL_OPERATOR::EXPLICIT_GPU:
        case FETIConfiguration::DUAL_OPERATOR::IMPLICIT_GPU:
        case FETIConfiguration::DUAL_OPERATOR::EXPLICIT_SC:
            eslog::error("not implemented dual operator\n");
            break;
        }
        break;
    }
    return dual;
}

template<typename T>
void DualOperator<T>::getInitVector(Vector_Dual<T> &v)
{
    double norm = 1 / std::sqrt(Dual_Map::total);
    math::set(v, T{0});
    for (size_t i = 0; i < Dual_Map::local_intervals.size(); ++i) {
        for (int j = Dual_Map::local_intervals[i].start; j < Dual_Map::local_intervals[i].end; ++j) {
            int mul = (Dual_Map::local_intervals[i].offset + j) % 2;
            v.vals[j] = (1 - 2 * mul) * norm;
        }
    }
    v.synchronize();
}

template<typename T>
void DualOperator<T>::reduceInfo(std::vector<DirectSparseSolver<T> > &KSolver, DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max)
{
    sum.nnzA = sum.nnzL = /*sum.memoryL =*/ sum.rows = sum.dualA = sum.surfaceA = 0;
    min.nnzA = min.nnzL = /*min.memoryL =*/ min.rows = min.dualA = min.surfaceA = INT32_MAX;
    max.nnzA = max.nnzL = /*max.memoryL =*/ max.rows = max.dualA = max.surfaceA = 0;
    for (size_t di = 0; di < feti.K.size(); ++di) {
        size_t dualA = feti.B1[di].nrows;
        size_t surfaceA = KSolver[di].getMatrixSize() - *std::min_element(feti.B1[di].cols, feti.B1[di].cols + feti.B1[di].nnz);
        min.rows = std::min<size_t>(min.rows, KSolver[di].getMatrixSize());
        min.nnzA = std::min<size_t>(min.nnzA, KSolver[di].getMatrixNnz());
        min.nnzL = std::min<size_t>(min.nnzL, KSolver[di].getFactorNnz());
        // min.memoryL = std::min(min.memoryL, KSolver[di].memoryL);
        min.dualA = std::min(min.dualA, dualA);
        min.surfaceA = std::min(min.surfaceA, surfaceA);
        max.rows = std::max<size_t>(max.rows, KSolver[di].getMatrixSize());
        max.nnzA = std::max<size_t>(max.nnzA, KSolver[di].getMatrixNnz());
        max.nnzL = std::max<size_t>(max.nnzL, KSolver[di].getFactorNnz());
        // max.memoryL = std::max(max.memoryL, KSolver[di].memoryL);
        max.dualA = std::max(max.dualA, dualA);
        max.surfaceA = std::max(max.surfaceA, surfaceA);
        sum.rows += KSolver[di].getMatrixSize();
        sum.nnzA += KSolver[di].getMatrixNnz();
        sum.nnzL += KSolver[di].getFactorNnz();
        // sum.memoryL += KSolver[di].memoryL;
        sum.dualA += dualA;
        sum.surfaceA += surfaceA;
    }

    Communication::allReduce(&min, nullptr, sizeof(DualOperatorInfo) / sizeof(size_t), MPITools::getType<size_t>().mpitype, MPI_MIN);
    Communication::allReduce(&max, nullptr, sizeof(DualOperatorInfo) / sizeof(size_t), MPITools::getType<size_t>().mpitype, MPI_MAX);
    Communication::allReduce(&sum, nullptr, sizeof(DualOperatorInfo) / sizeof(size_t), MPITools::getType<size_t>().mpitype, MPI_SUM);
}

template<typename T>
void DualOperator<T>::printInfo(std::vector<DirectSparseSolver<T> > &KSolver, DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max)
{
    eslog::info("      =  DOMAINS TOTAL                                                       %9d  = \n", feti.sinfo.domains);
    eslog::info("      =  DUAL SIZE                                                           %9d  = \n", Dual_Map::total);
    eslog::info("      =  B1 ROWS                                        %8.0f <%8d - %8d>  = \n", (double)sum.dualA / feti.sinfo.domains, min.dualA, max.dualA);
    eslog::info("      =  K+ SURFACE                                     %8.0f <%8d - %8d>  = \n", (double)sum.surfaceA / feti.sinfo.domains, min.surfaceA, max.surfaceA);
    eslog::info("      =  K+ ROWS                                        %8.0f <%8d - %8d>  = \n", (double)sum.rows / feti.sinfo.domains, min.rows, max.rows);
    eslog::info("      =  K+ NNZ                                         %8.0f <%8d - %8d>  = \n", (double)sum.nnzA / feti.sinfo.domains, min.nnzA, max.nnzA);
    eslog::info("      =  K+ FACTORS NNZ                                 %8.0f <%8d - %8d>  = \n", (double)sum.nnzL / feti.sinfo.domains, min.nnzL, max.nnzL);

// eslog::info(" =   K+ SOLVER MEMORY [MB]                                    %8.2f <%8.2f - %8.2f> = \n", (double)sum.memoryL / feti.sinfo.domains / 1024. / 1024., min.memoryL / 1024. / 1024., max.memoryL / 1024. / 1024.);
//    if (sparsity != DirectSparseSolver<T>::VectorSparsity::DENSE) {
//        eslog::info(" =   K+ FACTORIZATION                                                        RESPECT SURFACE = \n");
//    }
    if (feti.configuration.exhaustive_info) {
        T min = T{1e9}, max = T{0}, sum = T{0}, normK = T{0};
        for (size_t di = 0; di < feti.K.size(); ++di) {
            SpBLAS<Matrix_CSR, T> spblas(feti.K[di]);
            for (int r = 0; r < feti.K[di].nrows; ++r) {
                normK += feti.K[di].vals[feti.K[di].rows[r - Indexing::CSR]] * feti.K[di].vals[feti.K[di].rows[r - Indexing::CSR]];
            }
            normK = std::sqrt(normK);
            Vector_Dense<T> R, KR; KR.resize(feti.R1[di].ncols);
            T norm = T{0};
            for (int r = 0; r < feti.R1[di].nrows; ++r) {
                R.size = feti.R1[di].ncols; R.vals = feti.R1[di].vals + feti.R1[di].ncols * r;
                spblas.apply(KR, T{1}, T{0}, R);
                norm += math::norm(KR);
            }
            norm /= normK;
            min = std::min(min, norm);
            max = std::max(max, norm);
            sum += norm;
        }
        Communication::allReduce(&min, nullptr, 1, MPITools::getType<T>().mpitype, MPI_MIN);
        Communication::allReduce(&max, nullptr, 1, MPITools::getType<T>().mpitype, MPI_MAX);
        Communication::allReduce(&sum, nullptr, 1, MPITools::getType<T>().mpitype, MPI_SUM);
        eslog::info("      =  NORM(K * R) / NORM(DIAG(K))                    %.2e <%.2e - %.2e>  = \n", (double)sum / feti.sinfo.domains, min, max);

        // power method to Eigen values
        // B * Bt = eye
        // pseudo inverse
    }
}

template<typename T>
void DualOperator<T>::estimateMaxEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations)
{
    DualOperator<T> *F = feti.dualOperator;

    // {-1, 1} / norma
    Vector_Dual<T> y, v;
    getInitVector(y);

    lambda = std::sqrt(y.dot());
    double err = std::numeric_limits<T>::max();
    for (iterations = 0; iterations <= maxIterations && epsilon < err; ++iterations) {
        math::copy(v, y);
        math::scale(1 / lambda, v);
        F->apply(v, y);
        double _lambda = lambda;
        lambda = std::sqrt(y.dot());
        err = std::fabs(lambda - _lambda) / std::fabs(lambda);
    }
}

template<typename T>
void DualOperator<T>::estimateMaxProjectedEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations, double rho, double normPFP)
{
    DualOperator<T> *F = feti.dualOperator;
    Projector<double> *P = feti.projector;

    // {-1, 1} / norma
    Vector_Dual<T> v, y, x, z;
    getInitVector(y);
    lambda = std::sqrt(y.dot());
    double err = std::numeric_limits<T>::max();
    for (iterations = 1; iterations <= maxIterations && epsilon < err; ++iterations) {
        math::copy(x, y);
        math::scale(1 / lambda, x);

        // y = H(x)
        P->apply(x, v);
        F->apply(v, z);
        P->apply(z, y);
        double _lambda = lambda;
        lambda = std::sqrt(y.dot());
        err = std::fabs(lambda - _lambda) / std::fabs(lambda);
    }
}

template class DualOperator<double>;

}
