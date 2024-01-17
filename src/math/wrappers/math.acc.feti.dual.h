
#ifndef SRC_MATH_WRAPPERS_MATH_ACC_FETI_DUAL_H_
#define SRC_MATH_WRAPPERS_MATH_ACC_FETI_DUAL_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"
#include "feti/common/vector_dual.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename T>
struct Acc_FETI_Dual_Operator;

template <template <typename, typename> class Matrix, typename T, typename I = int>
struct AccFETIDualOperator {

	AccFETIDualOperator(int rank);
	~AccFETIDualOperator();

	void set(const std::vector<Matrix<T, I> > &K, const std::vector<Matrix<T, I> > &B);
	void update(const std::vector<Matrix<T, I> > &K);
	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y, const std::vector<std::vector<int> > & D2C);

private:
	int rank;
	Acc_FETI_Dual_Operator<T> *_acc;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_ACC_FETI_DUAL_H_ */
