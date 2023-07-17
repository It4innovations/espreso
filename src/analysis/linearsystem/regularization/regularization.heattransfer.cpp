
#include "regularization.heattransfer.h"

namespace espreso {

template struct RegularizationHeatTransfer<double>;

template <typename T>
void RegularizationHeatTransfer<T>::setAnalytic()
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		const Matrix_CSR<T> &K = feti.K.domains[d];
		Matrix_Dense<T> &R = feti.regularization.R1.domains[d];
		Matrix_CSR<T> &RegMat = feti.regularization.RegMat.domains[d];

		R.resize(K.nrows, 1);
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
void RegularizationHeatTransfer<T>::updateAnalytic()
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		const Matrix_CSR<T> &K = feti.K.domains[d];
		Matrix_Dense<T> &R = feti.regularization.R1.domains[d];
		Matrix_CSR<T> &RegMat = feti.regularization.RegMat.domains[d];

		RegMat.vals[0] = math::getDiagonalMax(K);
		math::set(R, 1.0 / std::sqrt(K.nrows));
	}

}

}
