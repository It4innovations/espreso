
#ifndef SRC_MATH_WRAPPERS_MATH_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SOLVER_H_

#include "math/primitives/matrix_dense.h"
#include "math/primitives/vector_dense.h"

#include <vector>

namespace espreso {

template <typename T, typename I = int>
struct DenseSolver {
    DenseSolver();
    DenseSolver(const Matrix_Dense<T, I> &a);
    DenseSolver(Matrix_Dense<T, I> &&a);
    DenseSolver(const DenseSolver &other) = delete;
    DenseSolver(DenseSolver &&other) = delete;
    DenseSolver & operator=(const DenseSolver &other) = delete;
    DenseSolver & operator=(DenseSolver &&other) = delete;
    ~DenseSolver();

    void commit(const Matrix_Dense<T, I> &a);
    void commit(Matrix_Dense<T, I> &&a);

    void factorization();

    void solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution);
    void solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution);

protected:
    Matrix_Dense<T, I> a;
    std::vector<I> ipiv;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_SOLVER_H_ */
