
#include "regularization.empty.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void RegularizationEmpty<T>::set(FETI<T> &feti)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		feti.R1[d].resize(0, feti.K[d].nrows);
		feti.R1[d].type = Matrix_Type::REAL_NONSYMMETRIC;
		feti.R1[d].shape = Matrix_Shape::FULL;

		feti.RegMat[d].resize(feti.K[d].nrows, feti.K[d].ncols, 0);
		feti.RegMat[d].type = feti.K[d].type;
		feti.RegMat[d].shape = feti.K[d].shape;
		for (esint r = 0; r <= feti.RegMat[d].nrows; ++r) {
			feti.RegMat[d].rows[r] = 0;
		}
	}
}

template <typename T>
void RegularizationEmpty<T>::update(FETI<T> &feti)
{

}

template struct RegularizationEmpty<double>;

}
