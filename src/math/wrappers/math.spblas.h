
#ifndef SRC_MATH_WRAPPERS_MATH_SPBLAS_H_
#define SRC_MATH_WRAPPERS_MATH_SPBLAS_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"
#include "math/primitives/matrix_ijv.h"

namespace espreso {

struct Matrix_SpBLAS_External_Representation;

template <template <typename, typename, typename> class Matrix, typename T, typename I = int>
struct SpBLAS {
    using MatrixType = Matrix<T, I, cpu_allocator>;

    SpBLAS();
    SpBLAS(MatrixType &a);
    ~SpBLAS();

    void insert(MatrixType &a);

    // y = alpha * A * x + beta * y
    void apply(Vector_Dense<T, I> &y, const T &alpha, const T &beta, const Vector_Dense<T, I> &x);

    MatrixType *matrix;

private:
    Matrix_SpBLAS_External_Representation *_spblas;
};

namespace math {
namespace spblas {

// input = output[start_row:end_row, start_col:end_col]. Start inclusive, end exclusive
template <typename T, typename I>
void submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, bool trans = false, bool conj = false, bool output_force_full = false);

// input = output[start_row:end_row, start_col:end_col]. Start inclusive, end exclusive
template <typename T, typename I>
void submatrix(const Matrix_CSR<T, I> &input, Matrix_CSR<T, I>   &output, I start_row, I end_row, I start_col, I end_col, bool trans = false, bool conj = false, bool output_force_full = false);

// TODO: do it in SpBLAS
// TODO: indexing
// y += alpha * A * x
// y += alpha * A * x
template <typename T, typename I>
void apply(Vector_Dense<T, I> &y, const T &alpha, Matrix_CSR<T, I> &A, const Vector_Dense<T> &x)
{
    if (A.rows[0] != 0) eslog::error("implement math::spblas::apply with non-zero based indexing.\n");
    for (int r = 0; r < A.nrows; ++r) {
        for (int c = A.rows[r]; c < A.rows[r + 1]; ++c) {
            y.vals[r] += alpha * A.vals[c] * x.vals[A.cols[c]];
        }
    }
}

template <typename T, typename I>
void applyT(Vector_Dense<T, I> &y, const T &alpha, Matrix_CSR<T, I> &A, const Vector_Dense<T> &x)
{
    if (A.rows[0] != 0) eslog::error("implement math::spblas::apply with non-zero based indexing.\n");
    for (int r = 0; r < A.nrows; ++r) {
        for (int c = A.rows[r]; c < A.rows[r + 1]; ++c) {
            y.vals[A.cols[c]] += alpha * A.vals[c] * x.vals[r];
        }
    }
}

template <typename T, typename I>
void apply(Vector_Dense<T, I> &y, const T &alpha, Matrix_CSR<T, I> &A, const int* D2C, const Vector_Dense<T> &x)
{
    if (A.rows[0] != 0) eslog::error("implement math::spblas::apply with non-zero based indexing.\n");
    for (int r = 0; r < A.nrows; ++r) {
        for (int c = A.rows[r]; c < A.rows[r + 1]; ++c) {
            y.vals[D2C[r]] += alpha * A.vals[c] * x.vals[A.cols[c]];
        }
    }
}

template <typename T, typename I>
void applyT(Vector_Dense<T, I> &y, const T &alpha, Matrix_CSR<T, I> &A, const int* D2C, const Vector_Dense<T> &x)
{
    if (A.rows[0] != 0) eslog::error("implement math::spblas::apply with non-zero based indexing.\n");
    for (int r = 0; r < A.nrows; ++r) {
        for (int c = A.rows[r]; c < A.rows[r + 1]; ++c) {
            y.vals[A.cols[c]] += alpha * A.vals[c] * x.vals[D2C[r]];
        }
    }
}

}
}

}

#include "math.spblas.hpp"

#endif /* SRC_MATH_WRAPPERS_MATH_SPBLAS_H_ */
