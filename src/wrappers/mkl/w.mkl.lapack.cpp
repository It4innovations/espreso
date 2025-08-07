
#include "math/math.h"
#include "esinfo/eslog.h"

#ifdef HAVE_MKL
#ifdef ESPRESO_USE_WRAPPER_LAPACK_MKL

#include "mkl_lapacke.h"

namespace espreso {
namespace math {
namespace lapack {

template <>
void solve_sym_upper(Matrix_Dense<double> &A, Matrix_Dense<double> &rhs)
{
    int *ipiv = new int[A.nrows];
    LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', A.nrows, A.vals, A.ncols, ipiv);
    LAPACKE_dsytrs(LAPACK_ROW_MAJOR, 'U', A.nrows, rhs.ncols, A.vals, A.ncols, ipiv, rhs.vals, rhs.ncols);
    delete[] ipiv;
}

template <>
void solve_general(Matrix_Dense<double> &A, Matrix_Dense<double> &rhs)
{
    int *ipiv = new int[A.nrows];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, A.nrows, rhs.ncols, A.vals, A.ncols, ipiv, rhs.vals, rhs.ncols);
    delete[] ipiv;
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values)
{
    int ret = 0;
    switch (A.shape) {
    case Matrix_Shape::FULL:  ret = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', A.nrows, A.vals, A.ncols, values.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'N', 'U', A.nrows, A.vals, values.vals, nullptr, A.ncols); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'N', 'L', A.nrows, A.vals, values.vals, nullptr, A.ncols); break;
    }
    if (ret) eslog::error("eigen values do not converge\n");
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, Matrix_Dense<double> &vectors)
{
    int ret = 0;
    // LAPACK_COL_MAJOR and swap 'L' and 'U' to get vectors in rows
    switch (A.shape) {
    case Matrix_Shape::FULL: math::copy(vectors, A); ret = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', A.nrows, vectors.vals, A.ncols, values.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'V', 'U', A.nrows, A.vals, values.vals, vectors.vals, A.ncols); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'V', 'L', A.nrows, A.vals, values.vals, vectors.vals, A.ncols); break;
    }
    if (ret) eslog::error("eigen values do not converge\n");
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, int begin, int end)
{
    int m; Vector_Dense<int, int> fail; fail.resize(2 * A.ncols);
    int ret = 0;
    switch (A.shape) {
    case Matrix_Shape::FULL:  ret = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'N', 'I', 'U', A.nrows, A.vals, A.ncols, 0, 0, begin, end, LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'N', 'I', 'U', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'N', 'I', 'L', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, nullptr, A.ncols, fail.vals); break;
    }
    if (ret) eslog::error("eigen values do not converge\n");
}

template <>
void get_eig_sym(Matrix_Dense<double> &A, Vector_Dense<double> &values, Matrix_Dense<double> &vectors, int begin, int end)
{
    int m; Vector_Dense<int, int> fail; fail.resize(2 * A.ncols);
    int ret = 0;
    switch (A.shape) {
    case Matrix_Shape::FULL:  ret = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', A.nrows, A.vals, A.ncols, 0, 0, begin, end, LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, vectors.ncols, fail.vals); break;
    case Matrix_Shape::UPPER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'V', 'I', 'U', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, A.ncols, fail.vals); break;
    case Matrix_Shape::LOWER: ret = LAPACKE_dspevx(LAPACK_ROW_MAJOR, 'V', 'I', 'L', A.nrows, A.vals, 0, 0, begin, end, 2 * LAPACKE_dlamch('S'), &m, values.vals, vectors.vals, A.ncols, fail.vals); break;
    }
    if (ret) eslog::error("eigen values do not converge\n");
}

template <>
void get_svd(Matrix_Dense<double> &A, Vector_Dense<double> &s, Matrix_Dense<double> &U, Matrix_Dense<double> &V)
{
    int ret = 0;
    s.resize(std::min(A.nrows, A.ncols)); math::set(s, 0.0);
    U.resize(A.nrows, A.nrows); math::set(U, 0.0);
    V.resize(A.ncols, A.ncols); math::set(V, 0.0);
    Vector_Dense<double> superb; superb.resize(s.size);
    ret = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', A.nrows, A.ncols, A.vals, A.ncols, s.vals, U.vals, U.ncols, V.vals, V.nrows, superb.vals);
    if (ret) eslog::error("error in 'get_svd'\n");
}

template <>
void submatrix(const Matrix_Dense<double, int> &input, Matrix_Dense<double, int> &output, int start_row, int end_row, int start_col, int end_col)
{
    if (input.shape != Matrix_Shape::FULL) eslog::error("Cannot copy triangular matrix.\n");

    output.resize(end_row - start_row, end_col - start_col);
    LAPACKE_dlacpy(LAPACK_ROW_MAJOR, ' ', end_row - start_row, end_col - start_col, input.vals + start_row * input.ncols + start_col, input.ncols, output.vals, output.ncols);
}

//void MATH::upDense3x3EigenValues(double *mVals, double *eigenValues)
//{
//#ifdef HAVE_MKL
//    double T[6] = { mVals[0], mVals[1], mVals[2], mVals[3], mVals[4], mVals[5] };
//    LAPACKE_dsterf(3, T, T + 3);
//    eigenValues[0] = T[2];
//    eigenValues[1] = T[1];
//    eigenValues[2] = T[0];
//#endif
//}
//
//void MATH::DenseMinGeneralizedEigenVectors(esint msize, double *A, double *B, esint n, double *lambdas, double *vectors)
//{
//#ifdef HAVE_MKL
//    esint nn;
//    esint *fail = new esint[msize];
//    double *tmpvec = new double[msize * msize];
//    LAPACKE_dsygvx(LAPACK_ROW_MAJOR, 1, 'V', 'I', 'U', msize, A, msize, B, msize, 0, 0, 1, n, 0, &nn, lambdas, tmpvec, msize, fail);
//    for (esint r = 0; r < msize; ++r) {
//        for (esint c = 0; c < n; ++c) {
//            vectors[r * n + c] = tmpvec[r * msize + c];
//        }
//    }
//    delete[] fail;
//    delete[] tmpvec;
//#endif
//}
//
//void MATH::upDense3x3EigenValuesEigenVectors(double *A, double *W, double *Z)
//{
//#ifdef HAVE_MKL
//    for(int i = 0; i < 9; ++i){
//        Z[i] = A[i];
//    }
//    esint info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', 3, Z, 3, W );
//
//    if(info > 0){
//        eslog::globalerror("ESPRESO internal error: upDense3x3EigenValuesEigenVectors failed to compued eigen-values or eigen-vectors.\n");
//    }
//#endif
//}
//
//void MATH::DenseMatDenseMatRowMajorSystemSolve(int nra, int ncb, double *a, double *b)
//{
//    // B = A\B
//#ifdef HAVE_MKL
//    //lapack_int LAPACKE_dgesv (int matrix_layout, lapack_int n, lapack_int nrhs, double * a, lapack_int lda, lapack_int * ipiv, double * b, lapack_int ldb)
//    esint ipiv[nra];
//    LAPACKE_dgesv(LAPACK_ROW_MAJOR, nra, ncb, a, nra, ipiv, b, ncb);
//#endif
//}

}
}
}

#endif
#endif
