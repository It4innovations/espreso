
#ifndef SRC_MATH_WRAPPERS_MATH_ACC_FETI_DUAL_H_
#define SRC_MATH_WRAPPERS_MATH_ACC_FETI_DUAL_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"
#include "math/feti/vector_dual.h"

#include <cstddef>
#include <vector>

namespace espreso {

struct Acc_FETI_Dual_Operator;

template <typename T, template <typename> class Matrix>
struct AccFETIDualOperator {

	AccFETIDualOperator();
	~AccFETIDualOperator();

	void set(const std::vector<Matrix_CSR<T> > &K, const std::vector<Matrix_CSR<T> > &B);
	void update(const std::vector<Matrix_CSR<T> > &K);
	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

private:
	Acc_FETI_Dual_Operator *_acc;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_ACC_FETI_DUAL_H_ */
