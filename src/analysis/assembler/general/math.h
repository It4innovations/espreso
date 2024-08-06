
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_

#include "wrappers/simd/simd.h"

namespace espreso {

template <size_t rows, size_t cols>
void set(SIMD M[rows * cols], const SIMD &v)
{
    for (size_t i = 0; i < rows * cols; ++i) {
        M[i] = v;
    }
}

template <size_t rows, size_t common, size_t cols>
void multAB(SIMD C[rows * cols], const SIMD A[rows * common], const SIMD B[common * cols], const SIMD &scale = load1(1.))
{
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            SIMD res;
            for (size_t k = 0; k < common; ++k) {
                res = res + A[i * common + k] * B[k * cols + j];
            }
            C[i * cols + j] = scale * res;
        }
    }
}

template <size_t rows, size_t common, size_t cols>
void multAtB(SIMD C[rows * cols], const SIMD A[rows * common], const SIMD B[common * cols], const SIMD &scale = load1(1.))
{
    for (size_t k = 0; k < common; ++k) {
        for (size_t j = 0; j < cols; ++j) {
            SIMD res;
            for (size_t i = 0; i < rows; ++i) {
                res = res + A[i * common + k] * B[i * cols + j];
            }
            C[k * cols + j] = C[k * cols + j] + scale * res;
        }
    }
}

template <size_t rows, size_t common, size_t cols>
void multABt(SIMD C[rows * cols], const SIMD A[rows * common], const SIMD B[common * cols], const SIMD &scale = load1(1.))
{
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            SIMD res;
            for (size_t k = 0; k < common; ++k) {
                res = res + A[i * common + k] * B[j * common + k];
            }
            C[i * cols + j] = C[i * cols + j] + scale * res;
        }
    }
}

template <size_t rows, size_t common, size_t cols>
void multAtBt(SIMD C[rows * cols], const SIMD A[rows * common], const SIMD B[common * cols], const SIMD &scale = load1(1.))
{
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            SIMD res;
            for (size_t k = 0; k < common; ++k) {
                res = res + A[k * rows + i] * B[j * common + k];
            }
            C[i * cols + j] = C[i * cols + j] + scale * res;
        }
    }
}

template <size_t rows, size_t cols>
void multAtBA(SIMD C[cols * cols], const SIMD A[rows * cols], const SIMD B[rows * rows], const SIMD &scale = load1(1.))
{
    for (size_t n = 0; n < cols; ++n) {
        SIMD AtB[rows];
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < rows; ++j) {
                AtB[j] = AtB[j] + A[i * cols + n] * B[i * rows + j];
            }
        }
        SIMD nn; for (size_t k = 0; k < rows; ++k) { nn = nn + AtB[k] * A[k * cols + n]; }
        C[n * cols + n] = C[n * cols + n] + scale * nn;
        for (size_t m = n + 1; m < cols; ++m) {
            SIMD nm; for (size_t k = 0; k < rows; ++k) { nm = nm + AtB[k] * A[k * cols + m]; }
            C[n * cols + m] = C[n * cols + m] + scale * nm;
            C[m * cols + n] = C[m * cols + n] + scale * nm;
        }
    }
}

template <size_t rows, size_t cols>
void multABAt(SIMD C[cols * cols], const SIMD A[rows * cols], const SIMD B[rows * rows], const SIMD &scale = load1(1.))
{
    for (size_t n = 0; n < rows; ++n) {
        SIMD AB[cols];
        for (size_t i = 0; i < cols; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                AB[j] = AB[j] + A[n * cols + i] * B[i * cols + j];
            }
        }
        SIMD nn; for (size_t k = 0; k < cols; ++k) { nn = nn + AB[k] * A[n * cols + k]; }
        C[n * rows + n] = C[n * rows + n] + scale * nn;
        for (size_t m = n + 1; m < rows; ++m) {
            SIMD nm; for (size_t k = 0; k < cols; ++k) { nm = nm + AB[k] * A[m * cols + k]; }
            C[n * rows + m] = C[n * rows + m] + scale * nm;
            C[m * rows + n] = C[m * rows + n] + scale * nm;
        }
    }
}


inline SIMD determinant(const SIMD J[9])
{
    return
            + J[0] * J[4] * J[8] + J[1] * J[5] * J[6] + J[2] * J[3] * J[7]
            - J[2] * J[4] * J[6] - J[1] * J[3] * J[8] - J[0] * J[5] * J[7];
}

void eigSym(const SIMD A[9], SIMD eVal[3]);
void eigSym(const SIMD A[9], SIMD eVal[3], SIMD eVec[9]);

inline void inv(const SIMD A[9], SIMD &det, SIMD invA[9])
{
    det = determinant(A);
    SIMD scale = ones() / det;
    invA[0] = scale * ( A[8] * A[4] - A[7] * A[5]);
    invA[1] = scale * (-A[8] * A[1] + A[7] * A[2]);
    invA[2] = scale * ( A[5] * A[1] - A[4] * A[2]);
    invA[3] = scale * (-A[8] * A[3] + A[6] * A[5]);
    invA[4] = scale * ( A[8] * A[0] - A[6] * A[2]);
    invA[5] = scale * (-A[5] * A[0] + A[3] * A[2]);
    invA[6] = scale * ( A[7] * A[3] - A[6] * A[4]);
    invA[7] = scale * (-A[7] * A[0] + A[6] * A[1]);
    invA[8] = scale * ( A[4] * A[0] - A[3] * A[1]);
}

inline void inv(const SIMD A[9], SIMD invA[9])
{
    SIMD det; inv(A, det, invA);
}


}

#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_ */
