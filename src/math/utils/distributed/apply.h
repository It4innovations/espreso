
#ifndef SRC_MATH2_UTILS_DISTRIBUTED_APPLY_H_
#define SRC_MATH2_UTILS_DISTRIBUTED_APPLY_H_

#include "math/primitives/matrix_csr.h"
#include "math/wrappers/math.spblas.h"

#include <vector>

namespace espreso {

template <typename T> class Vector_Dense;
template <template<typename> typename Matrix, typename T> class Vector_Distributed;
template <template<typename> typename Matrix, typename T> class Matrix_Distributed;

template <template<typename> typename Struct, typename T> struct Data_Apply { };

template <typename T>
struct Data_Apply<Matrix_CSR, T> {
	Matrix_CSR<T> m;
	Vector_Dense<T> v;
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset, sOffset;
	std::vector<int> neighbors;
	std::vector<esint> nDOF;
	esint offset;

	SpBLAS<T, Matrix_CSR> spblas;

	Data_Apply<Matrix_CSR, T>& operator=(const Data_Apply<Matrix_CSR, T> &other)
	{
		m._Matrix_CSR<T>::operator=(other.m);
		return *this;
	}

	void init(Matrix_Distributed<Matrix_CSR, T> &m);
	void commit(Matrix_Distributed<Matrix_CSR, T> &m);
	void apply(Vector_Distributed<Vector_Dense, T> *y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> *x);
};

}

#endif /* SRC_MATH2_UTILS_DISTRIBUTED_APPLY_H_ */
