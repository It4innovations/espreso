
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

        if (feti.configuration.projector == FETIConfiguration::PROJECTOR::CONJUGATE) {
            feti.RegMat[d].resize(feti.K[d].nrows, feti.K[d].ncols, 0);
            feti.RegMat[d].type = feti.K[d].type;
            feti.RegMat[d].shape = feti.K[d].shape;
            for (esint r = 0; r <= feti.RegMat[d].nrows; ++r) {
                feti.RegMat[d].rows[r] = 0;
            }
        } else {
            RegMat.resize(K.nrows, K.ncols, 1);
            RegMat.type = K.type;
            RegMat.shape = K.shape;

            RegMat.rows[0] = RegMat.cols[0] = Indexing::CSR;
            std::fill(RegMat.rows + 1, RegMat.rows + RegMat.nrows + 1, RegMat.rows[0] + 1);
        }
    }
}

template <typename T>
void RegularizationHeatTransfer<T>::update(FETI<T> &feti)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        math::set(feti.R1[d], T{1});
        if (feti.configuration.projector != FETIConfiguration::PROJECTOR::CONJUGATE) {
            feti.RegMat[d].vals[0] = 0;
            for (esint r = 0; r < feti.K[d].nrows; ++r) {
                feti.RegMat[d].vals[0] = std::max(feti.RegMat[d].vals[0], feti.K[d].vals[feti.K[d].rows[r] - Indexing::CSR]);
            }
        }
    }
}

template struct RegularizationHeatTransfer<double>;

}
