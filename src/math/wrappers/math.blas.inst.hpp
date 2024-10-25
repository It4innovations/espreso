
namespace espreso {
namespace math {
namespace blas {

    #define INSTANTIATE_T_I(T,I) \
    template void apply(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x); \
    template void applyT(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x); \
    template void apply_hermitian(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, char uplo, const T &beta, const Vector_Dense<T, I> &x); \
    template void AAt(const Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &AAt, bool trans); \
    template void multiply(T alpha, const Matrix_Dense<T, I> &A, const Matrix_Dense<T, I> &B, T beta, Matrix_Dense<T, I> &C, bool transA, bool transB); \
    template void multiply(T alpha, const Matrix_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, Vector_Dense<T, I> &C, bool transA); \
    template void multiply(T alpha, const Vector_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, T &out);

        #define INSTANTIATE_T(T) \
        INSTANTIATE_T_I(T, int32_t) \
        /* INSTANTIATE_T_I(T, int64_t) */ \
        template void copy(const int size, T *x, const int incX, const T *y, const int incY); \
        template void scale(const int size, const T &alpha, T *x, const int incX); \
        template void add(const int size, T *x, const int incX, const T &alpha, const T *y, const int incY); \
        template T dot(const int size, const T *x, const int incX, const T *y, const int incY); \
        template utils::remove_complex_t<T> norm(const int size, const T *x, const int incX);

            INSTANTIATE_T(float)
            INSTANTIATE_T(double)
            /* INSTANTIATE_T(std::complex<float>) */
            INSTANTIATE_T(std::complex<double>)

        #undef INSTANTIATE_T
    #undef INSTANTIATE_T_I

}
}
}
