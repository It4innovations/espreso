
#ifndef SRC_ANALYSIS_BUILDER_DIRECT_APPLY_H_
#define SRC_ANALYSIS_BUILDER_DIRECT_APPLY_H_

#include "math/wrappers/math.spblas.h"

#include <vector>

namespace espreso {

template <typename T> class Matrix_Distributed;
template <template<typename, typename, typename> typename Vector, typename T> class Vector_Distributed;

template <typename T>
struct Matrix_CSR_Apply {
	Matrix_CSR<T> localM;
	Vector_Dense<T> localV;
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset, sOffset;
	std::vector<int> neighbors;
	std::vector<esint> nDOF;
	esint offset;

	SpBLAS<Matrix_CSR, T> spblas;

	void init(Matrix_Distributed<T> &m);

	void apply(Matrix_Distributed<T> &m, Vector_Distributed<Vector_Dense, T> *y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> *x);
};

}



#endif /* SRC_ANALYSIS_BUILDER_DIRECT_APPLY_H_ */
