
#include "regularization.heattransfer.h"
#include "math/math.h"
#include "esinfo/ecfinfo.h"

namespace espreso {

template <typename T>
void RegularizationHeatTransfer<T>::set(FETI<T> &feti)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        const Matrix_CSR<T> &K = feti.K[d];
        Matrix_Dense<T> &R = feti.R1[d];
        Matrix_CSR<T> &RegMat = feti.RegMat[d];

        R.resize(1, K.nrows);
        R.type = Matrix_Type::REAL_NONSYMMETRIC;
        R.shape = Matrix_Shape::FULL;

        RegMat.resize(K.nrows, K.ncols, 1);
        RegMat.type = K.type;
        RegMat.shape = K.shape;

        RegMat.rows[0] = RegMat.cols[0] = Indexing::CSR;
        std::fill(RegMat.rows + 1, RegMat.rows + RegMat.nrows + 1, RegMat.rows[0] + 1);
    }


}

template <typename T>
void RegularizationHeatTransfer<T>::update(FETI<T> &feti)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        const Matrix_CSR<T> &K = feti.K[d];
        Matrix_Dense<T> &R = feti.R1[d];
        Matrix_CSR<T> &RegMat = feti.RegMat[d];

        RegMat.vals[0] = 0;
        for (esint r = 0; r < K.nrows; ++r) {
            RegMat.vals[0] = std::max(RegMat.vals[0], K.vals[K.rows[r] - Indexing::CSR]);
        }
        math::set(R.ncols, R.vals, 1, 1.0 / std::sqrt(K.nrows));
    }
}

template struct RegularizationHeatTransfer<double>;

}
