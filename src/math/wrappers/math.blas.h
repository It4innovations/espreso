
#ifndef SRC_MATH_WRAPPERS_MATH_BLAS_H_
#define SRC_MATH_WRAPPERS_MATH_BLAS_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_ijv.h"
#include "basis/utilities/utils.h"
#include "math/primitives_new.h"
#include "math/primitives_new/allocator_new.h"

#include <complex>

namespace espreso {
namespace math {
namespace blas {

    enum struct herk_mode { AhA, AAh };

    // x = y
    template <typename T>
    void copy(const int size, T *x, const int incX, const T *y, const int incY);

    // x *= alpha
    template <typename T>
    void scale(const int size, const T &alpha, T *x, const int incX);

    // x += alpha * y
    template <typename T>
    void add(const int size, T *x, const int incX, const T &alpha, const T *y, const int incY);

    template <typename T>
    T dot(const int size, const T *x, const int incX, const T *y, const int incY);

    template <typename T>
    utils::remove_complex_t<T> norm(const int size, const T *x, const int incX);

    // y = alpha * A * x + beta * y
    template <typename T, typename I>
    void apply(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x);

    // y = alpha * At * x + beta * y
    template <typename T, typename I>
    void applyT(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x);

    template <typename T, typename I>
    void apply_hermitian(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, char uplo, const T &beta, const Vector_Dense<T, I> &x);

    template <typename T, typename I>
    void AAt(const Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &AAt, bool trans = false);

    // C = alpha * op(A) * op(B) + beta C
    template <typename T, typename I>
    void multiply(T alpha, const Matrix_Dense<T, I> &A, const Matrix_Dense<T, I> &B, T beta, Matrix_Dense<T, I> &C, bool transA = false, bool transB = false);

    template <typename T, typename I>
    void multiply(T alpha, const Matrix_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, Vector_Dense<T, I> &C, bool transA = false);

    template <typename T, typename I>
    void multiply(T alpha, const Vector_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, T &out);

    template<typename T>
    void trsm(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & X, T alpha = T{1});

    template<typename T>
    void gemm(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & B, MatrixDenseView_new<T> & C, T alpha = T{1}, T beta = T{0});

    template<typename T>
    void herk(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & C, herk_mode mode, utils::remove_complex_t<T> alpha, utils::remove_complex_t<T> beta);

    template<typename T, typename I>
    void matrix_conj(T * A, I nrows, I ncols, I ld, char order, char uplo)
    {
        static_assert(utils::is_complex<T>(), "only complex types supported");
        I size_primary = ((order == 'R') ? nrows : ncols);
        I size_secdary = ((order == 'R') ? ncols : nrows);
        utils::remove_complex_t<T> * A_real = reinterpret_cast<utils::remove_complex_t<T>*>(A);
        for(I ip = 0; ip < size_primary; ip++) {
            I start_secdary = 0;
            I end_secdary = size_secdary;
            if((uplo == 'U' && order == 'R') || (uplo == 'L' && order == 'C')) start_secdary = ip;
            if((uplo == 'U' && order == 'C') || (uplo == 'L' && order == 'R')) end_secdary = ip;
            size_t use_size_secdary = end_secdary - start_secdary;
            scale(use_size_secdary, utils::remove_complex_t<T>{-1}, A_real + 2 * ip * ld + 2 * start_secdary + 1, 2);
        }
    }

    template<typename T>
    void transpose(size_t src_nrows, size_t src_ncols, const T * src, size_t src_ld, T * dst, size_t dst_ld, char order, bool conj);

    template<typename T>
    void transpose_inplace(size_t size, T * matrix, size_t ld, char order, bool conj);

    template <typename T, typename I>
    void apply(VectorDenseView_new<T> &y, const T &alpha, const MatrixDenseView_new<T> &a, const T &beta, const VectorDenseView_new<T> &x)
    {
        Vector_Dense<T,I> y_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(y);
        Vector_Dense<T,I> x_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(x);

        if(a.order == 'C') {
            MatrixDenseView_new<T> At = a.get_transposed_reordered_view();
            Matrix_Dense<T,I> At_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(At);

            applyT<T,I>(y_old, alpha, At_old, beta, x_old);
        }
        if(a.order == 'R') {
            Matrix_Dense<T,I> A_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(a);

            apply<T,I>(y_old, alpha, A_old, beta, x_old);
        }
    }

    template <typename T, typename I>
    void apply_hermitian(VectorDenseView_new<T> &y, const T &alpha, const MatrixDenseView_new<T> &a, const T &beta, const VectorDenseView_new<T> &x)
    {
        if(!is_hermitian<T>(a.prop.symm)) eslog::error("matrix has to be hermitian\n");
        if(a.prop.uplo != 'U' && a.prop.uplo != 'L') eslog::error("unset uplo\n");

        if(a.order != 'R') {
            MatrixDenseView_new<T> At = a.get_transposed_reordered_view();
            if constexpr(utils::is_complex<T>()) {
                VectorDenseData_new<T> x_conj;
                x_conj.set(x.size, AllocatorCPU_new::get_singleton());
                x_conj.alloc();
                utils::transform_n(x.vals, x.size, x_conj.vals, [](const T & x){return std::conj(x);});
                apply_hermitian<T,I>(y, alpha, At, beta, x);
                utils::transform_n(y.vals, y.size, y.vals, [](const T & x){return std::conj(x);});
            }
            else {
                apply_hermitian<T,I>(y, alpha, At, beta, x);
            }
            return;
        }

        Vector_Dense<T,I> y_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(y);
        Vector_Dense<T,I> x_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(x);
        Matrix_Dense<T,I> A_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(a);
        char uplo = a.prop.uplo;

        apply_hermitian<T,I>(y_old, alpha, A_old, uplo, beta, x_old);
    }

    template<typename T>
    void hemm(size_t m, size_t n, T alpha, T * A, char order_A, size_t lda, char uplo_A, T * B, char order_B, size_t ldb, T beta, T * C, char order_C, size_t ldc);

    template<typename T>
    void hemm(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & B, MatrixDenseView_new<T> & C, T alpha = T{1}, T beta = T{0})
    {
        hemm(C.nrows, C.ncols, alpha, A.vals, A.order, A.ld, A.prop.uplo, B.vals, B.order, B.ld, beta, C.vals, C.order, C.ld);
    }

}
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_BLAS_H_ */
