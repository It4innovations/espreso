
#ifndef SRC_FETI_COMMON_MATRIX_DUAL_SPARSE_H_
#define SRC_FETI_COMMON_MATRIX_DUAL_SPARSE_H_

#include "matrix_dual.h"

namespace espreso {

template <typename T>
struct Matrix_Dual_Sparse: public Matrix_Dual<T> {

    Matrix_Dual_Sparse(const std::vector<int> &nonzeros, const std::vector<int> &distributed);

    virtual void synchronize() override;

    using Matrix_Dense<T>::nrows;
    using Matrix_Dense<T>::ncols;
    using Matrix_Dense<T>::vals;

    std::vector<int> nonzeros, distributed;
    std::vector<std::vector<int> > recvNonzeros;
};

}

#endif /* SRC_FETI_COMMON_MATRIX_DUAL_SPARSE_H_ */
