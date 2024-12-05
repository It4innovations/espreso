
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
void eigSym22Desc(const SIMD A[4], SIMD eVal[2]);

void eigSym33Desc(const SIMD A[9], SIMD eVal[3]);
void eigSym33Desc(const SIMD A[9], SIMD eVal[3], SIMD eVec[9]);
void eigSym33Asc(const SIMD A[9], SIMD eVal[3], SIMD eVec[9]);

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

inline void voigt3ToMatrix22(const SIMD voigt[3], SIMD A[4])
{
    A[0] = voigt[0]; A[1] = voigt[2];
    A[2] = voigt[2]; A[3] = voigt[1];
}

inline int voigt36ToMatrix3333(int i)
{
    switch (i) {
    case  0: return  0;
    case  1: return  4;
    case  2: return  8;
    case  3: return  1;
    case  4: return  5;
    case  5: return  2;
    case  6: return  4;
    case  7: return 40;
    case  8: return 44;
    case  9: return 37;
    case 10: return 41;
    case 11: return 38;
    case 12: return  8;
    case 13: return 44;
    case 14: return 80;
    case 15: return 73;
    case 16: return 77;
    case 17: return 74;
    case 18: return  1;
    case 19: return 37;
    case 20: return 73;
    case 21: return 10;
    case 22: return 14;
    case 23: return 11;
    case 24: return  5;
    case 25: return 41;
    case 26: return 77;
    case 27: return 14;
    case 28: return 50;
    case 29: return 47;
    case 30: return  2;
    case 31: return 38;
    case 32: return 74;
    case 33: return 11;
    case 34: return 47;
    case 35: return 20;
    }
    return -1;
}

inline int matrix33ToVoigh6(int i)
{
    switch (i) {
    case  0: return  0;
    case  1: return  3;
    case  2: return  5;
    case  3: return  3;
    case  4: return  1;
    case  5: return  4;
    case  6: return  5;
    case  7: return  4;
    case  8: return  2;
    }
    return -1;
}

inline int voigt6ToMatrix33(int i)
{
    switch (i) {
    case  0: return  0;
    case  1: return  4;
    case  2: return  8;
    case  3: return  1;
    case  4: return  5;
    case  5: return  2;
    }
    return -1;
}

inline int matrix3333ToVoigh36(int i)
{
    switch (i) {
    case  0: return  0;
    case  1: return  3;
    case  2: return  5;
    case  3: return  3;
    case  4: return  1;
    case  5: return  4;
    case  6: return  5;
    case  7: return  4;
    case  8: return  2;
    case  9: return  3;
    case 10: return 21;
    case 11: return 23;
    case 12: return 21;
    case 13: return  9;
    case 14: return 22;
    case 15: return 23;
    case 16: return 22;
    case 17: return 15;
    case 18: return  5;
    case 19: return 23;
    case 20: return 35;
    case 21: return 23;
    case 22: return 11;
    case 23: return 29;
    case 24: return 35;
    case 25: return 29;
    case 26: return 17;
    case 27: return  3;
    case 28: return 21;
    case 29: return 23;
    case 30: return 21;
    case 31: return  9;
    case 32: return 22;
    case 33: return 23;
    case 34: return 22;
    case 35: return 15;
    case 36: return  1;
    case 37: return  9;
    case 38: return 11;
    case 39: return  9;
    case 40: return  7;
    case 41: return 10;
    case 42: return 11;
    case 43: return 10;
    case 44: return  8;
    case 45: return  4;
    case 46: return 22;
    case 47: return 29;
    case 48: return 22;
    case 49: return 10;
    case 50: return 28;
    case 51: return 29;
    case 52: return 28;
    case 53: return 16;
    case 54: return  5;
    case 55: return 23;
    case 56: return 35;
    case 57: return 23;
    case 58: return 11;
    case 59: return 29;
    case 60: return 35;
    case 61: return 29;
    case 62: return 17;
    case 63: return  4;
    case 64: return 22;
    case 65: return 29;
    case 66: return 22;
    case 67: return 10;
    case 68: return 28;
    case 69: return 29;
    case 70: return 28;
    case 71: return 16;
    case 72: return  2;
    case 73: return 15;
    case 74: return 17;
    case 75: return 15;
    case 76: return  8;
    case 77: return 16;
    case 78: return 17;
    case 79: return 16;
    case 80: return 14;
    }
    return -1;
}

void print(size_t rows, size_t cols, const SIMD *A);


}

#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_ */
