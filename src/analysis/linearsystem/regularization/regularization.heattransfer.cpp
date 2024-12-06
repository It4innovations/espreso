
#include "regularization.h"
#include "math/math.h"
#include "esinfo/ecfinfo.h"

namespace espreso {

template <typename T>
static void setR1(Matrix_CSR<T> &K, Matrix_Dense<T> &R1)
{
    R1.resize(1, K.nrows);
    R1.type = Matrix_Type::REAL_NONSYMMETRIC;
    R1.shape = Matrix_Shape::FULL;
}

template <typename T>
static void updateR1(Matrix_CSR<T> &K, Matrix_Dense<T> &R1)
{
    math::set(R1, T{1});
}

template <typename T>
static void setRegMat(Matrix_CSR<T> &K, Matrix_CSR<T> &RegMat)
{
    RegMat.resize(K.nrows, K.ncols, 1);
    RegMat.type = K.type;
    RegMat.shape = K.shape;

    RegMat.rows[0] = RegMat.cols[0] = Indexing::CSR;
    std::fill(RegMat.rows + 1, RegMat.rows + RegMat.nrows + 1, RegMat.rows[0] + 1);
}

template <typename T>
static void updateRegMat(Matrix_CSR<T> &K, Matrix_CSR<T> &RegMat)
{
    RegMat.vals[0] = 0;
    for (esint r = 0; r < K.nrows; ++r) {
        RegMat.vals[0] = std::max(RegMat.vals[0], K.vals[K.rows[r] - Indexing::CSR]);
    }
}

template <typename T>
void Regularization<T>::set(FETI<T> &feti, HeatTransferLoadStepConfiguration &configuration)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        if (R1)     setR1    (feti.K[d], feti.R1[d]);
        if (regMat) setRegMat(feti.K[d], feti.RegMat[d]);
    }
}

template <typename T>
void Regularization<T>::update(FETI<T> &feti, HeatTransferLoadStepConfiguration &configuration)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        if (R1     && feti.updated.K) updateR1    (feti.K[d], feti.R1[d]);
        if (regMat && feti.updated.K) updateRegMat(feti.K[d], feti.RegMat[d]);
    }
    if (R1 && feti.updated.K) {
        orthonormalize(feti);
    }
}

template void Regularization<double>::set    (FETI<double> &feti, HeatTransferLoadStepConfiguration &configuration);
template void Regularization<double>::update (FETI<double> &feti, HeatTransferLoadStepConfiguration &configuration);

}
