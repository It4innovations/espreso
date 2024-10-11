
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
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            SIMD res;
            for (size_t k = 0; k < common; ++k) {
                res = res + A[r * common + k] * B[k * cols + c];
            }
            C[r * cols + c] = C[r * cols + c] + scale * res;
        }
    }
}

template <size_t rows, size_t common, size_t cols>
void multAtB(SIMD C[rows * cols], const SIMD A[common * rows], const SIMD B[common * cols], const SIMD &scale = load1(1.))
{
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            SIMD res;
            for (size_t k = 0; k < common; ++k) {
                res = res + A[k * rows + r] * B[k * cols + c];
            }
            C[r * cols + c] = C[r * cols + c] + scale * res;
        }
    }
}

template <size_t rows, size_t common, size_t cols>
void multABt(SIMD C[rows * cols], const SIMD A[rows * common], const SIMD B[cols * common], const SIMD &scale = load1(1.))
{
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            SIMD res;
            for (size_t k = 0; k < common; ++k) {
                res = res + A[r * common + k] * B[c * common + k];
            }
            C[r * cols + c] = C[r * cols + c] + scale * res;
        }
    }
}

template <size_t rows, size_t common, size_t cols>
void multAtBt(SIMD C[rows * cols], const SIMD A[rows * common], const SIMD B[common * cols], const SIMD &scale = load1(1.))
{
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            SIMD res;
            for (size_t k = 0; k < common; ++k) {
                res = res + A[k * rows + r] * B[c * common + k];
            }
            C[r * cols + c] = C[r * cols + c] + scale * res;
        }
    }
}

template <size_t rows, size_t cols>
void multAtBA(SIMD C[cols * cols], const SIMD A[rows * cols], const SIMD B[rows * rows], const SIMD &scale = load1(1.))
{
    for (size_t r = 0; r < cols; ++r) {
        SIMD AtB[rows];
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < rows; ++j) {
                AtB[i] = AtB[i] + A[j * cols + r] * B[j * rows + i];
            }
        }
        for (size_t c = 0; c < cols; ++c) {
            SIMD res;
            for (size_t k = 0; k < rows; ++k) {
                res = res + AtB[k] * A[k * cols + c];
            }
            C[r * cols + c] = C[r * cols + c] + scale * res;
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
        for (size_t m = 0; m < rows; ++m) {
            SIMD nm;
            for (size_t k = 0; k < cols; ++k) {
                nm = nm + AB[k] * A[m * cols + k];
            }
            C[n * rows + m] = C[n * rows + m] + scale * nm;
        }
    }
}


void eigSym22(const SIMD A[4], SIMD eVal[2]);
void eigSym22(const SIMD A[4], SIMD eVal[2], SIMD eVec[4]);

void eigSym33(const SIMD A[9], SIMD eVal[3]);
void eigSym33(const SIMD A[9], SIMD eVal[3], SIMD eVec[9]);

inline SIMD determinant22(const SIMD J[4])
{
    return J[0] * J[3] - J[1] * J[2];
}

inline void inv22(const SIMD A[4], SIMD &det, SIMD invA[4])
{
    det = determinant22(A);
    SIMD scale = ones() / det;
    invA[0] =  scale * A[3];
    invA[1] = -scale * A[1];
    invA[2] = -scale * A[2];
    invA[3] =  scale * A[0];
}

inline void inv22(const SIMD A[4], SIMD invA[4])
{
    SIMD det; inv22(A, det, invA);
}

inline SIMD determinant33(const SIMD J[9])
{
    return
            + J[0] * J[4] * J[8] + J[1] * J[5] * J[6] + J[2] * J[3] * J[7]
            - J[2] * J[4] * J[6] - J[1] * J[3] * J[8] - J[0] * J[5] * J[7];
}

inline void inv33(const SIMD A[9], SIMD &det, SIMD invA[9])
{
    det = determinant33(A);
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

inline void inv33(const SIMD A[9], SIMD invA[9])
{
    SIMD det; inv33(A, det, invA);
}

inline void matrix33ToVoigt6(const SIMD A[9], SIMD voigt[6])
{
    voigt[0] = A[0];
    voigt[1] = A[4];
    voigt[2] = A[8];
    voigt[3] = A[1];
    voigt[4] = A[5];
    voigt[5] = A[2];
}

inline void voigt6ToMatrix33(const SIMD voigt[6], SIMD A[9])
{
    A[0] = voigt[0]; A[1] = voigt[3]; A[2] = voigt[5];
    A[3] = voigt[3]; A[4] = voigt[1]; A[5] = voigt[4];
    A[6] = voigt[5]; A[7] = voigt[4]; A[8] = voigt[2];
}

void print(size_t rows, size_t cols, const SIMD *A);


}

#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_ */
