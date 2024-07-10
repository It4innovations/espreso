
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_

#include "wrappers/simd/simd.h"

namespace espreso {

template <size_t dim>
void multAB(SIMD C[dim * dim], const SIMD A[dim * dim], const SIMD B[dim * dim], const double &scaleAB)
{
    SIMD ab = load1(scaleAB);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            for (size_t k = 0; k < dim; ++k) {
                C[i * dim + j] = C[i * dim + j] + ab * A[i * dim + k] * B[k * dim + j];
            }
        }
    }
}

template <size_t dim>
void multABt(SIMD C[dim * dim], const SIMD A[dim * dim], const SIMD B[dim * dim], const double &scaleAB)
{
    SIMD ab = load1(scaleAB);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            for (size_t k = 0; k < dim; ++k) {
                C[i * dim + j] = C[i * dim + j] + ab * A[i * dim + k] * B[j * dim + k];
            }
        }
    }
}

template <size_t dim>
void multAtB(SIMD C[dim * dim], const SIMD A[dim * dim], const SIMD B[dim * dim], const double &scaleAB)
{
    SIMD ab = load1(scaleAB);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            for (size_t k = 0; k < dim; ++k) {
                C[i * dim + j] = C[i * dim + j] + ab * A[k * dim + i] * B[k * dim + j];
            }
        }
    }
}

template <size_t dim>
void multAtBt(SIMD C[dim * dim], const SIMD A[dim * dim], const SIMD B[dim * dim], const double &scaleAB)
{
    SIMD ab = load1(scaleAB);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            for (size_t k = 0; k < dim; ++k) {
                C[i * dim + j] = C[i * dim + j] + ab * A[k * dim + i] * B[j * dim + k];
            }
        }
    }
}

template <size_t rows, size_t cols>
void multMv(SIMD C[rows], const SIMD A[rows * cols], const SIMD B[cols], const SIMD &scale)
{
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            C[r] = C[r] + scale * A[r * cols + c] * B[c];
        }
    }
}

template <size_t rows, size_t cols>
void multMtv(SIMD C[cols], const SIMD A[rows * cols], const SIMD B[rows], const SIMD &scale)
{
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            C[c] = C[c] + scale * A[r * cols + c] * B[r];
        }
    }
}

template <size_t rows, size_t cols>
void multAtBA(SIMD C[cols * cols], const SIMD A[rows * cols], const SIMD B[rows * rows], const SIMD &scale)
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

}

#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_MATH_H_ */
