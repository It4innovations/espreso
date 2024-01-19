
#ifndef SRC_ANALYSIS_LINEARSYSTEM_MATRICES_MATRIX_DISTRIBUTED_APPLY_H_
#define SRC_ANALYSIS_LINEARSYSTEM_MATRICES_MATRIX_DISTRIBUTED_APPLY_H_

#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/vector_dense.h"
#include "math/wrappers/math.spblas.h"

#include <vector>

namespace espreso {

template <template<typename, typename, template<typename> typename> typename Vector, typename T> class Vector_Distributed;
template <template<typename, typename, template<typename> typename> typename Matrix, typename T> class Matrix_Distributed;

template <template<typename, typename, template<typename> typename> typename Struct, typename T, typename I = esint> struct Data_Apply { };

template <typename T, typename I>
struct Data_Apply<Matrix_CSR, T, I> {
	Matrix_CSR<T, I> m;
	Vector_Dense<T, I> v;
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<I> > rOffset, sOffset;
	std::vector<int> neighbors;
	std::vector<I> nDOF;
	I offset;

	SpBLAS<Matrix_CSR, T, I> spblas;

	Data_Apply<Matrix_CSR, T>& operator=(const Data_Apply<Matrix_CSR, T> &other)
	{
		m._Matrix_CSR<T, I>::operator=(other.m);
		return *this;
	}

	void init(Matrix_Distributed<Matrix_CSR, T> &m);
	void commit(Matrix_Distributed<Matrix_CSR, T> &m);
	void apply(Vector_Distributed<Vector_Dense, T> *y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> *x);
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_MATRICES_MATRIX_DISTRIBUTED_APPLY_H_ */
