
#include "math.h"

#include <cstdio>

namespace espreso {

void eigSym(const SIMD A[9], SIMD eVal[3])
{
    // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
    // p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2
    // if (p1 == 0)
    //    % A is diagonal.
    //    eig1 = A(1,1)
    //    eig2 = A(2,2)
    //    eig3 = A(3,3)
    // else
    //    q = trace(A)/3               % trace(A) is the sum of all diagonal values
    //    p2 = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1
    //    p = sqrt(p2 / 6)
    //    B = (1 / p) * (A - q * I)    % I is the identity matrix
    //    r = det(B) / 2
    //
    //    % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
    //    % but computation error can leave it slightly outside this range.
    //    if (r <= -1)
    //       phi = pi / 3
    //    elseif (r >= 1)
    //       phi = 0
    //    else
    //       phi = acos(r) / 3
    //    end
    //
    //    % the eigenvalues satisfy eig3 <= eig2 <= eig1
    //    eig1 = q + 2 * p * cos(phi)
    //    eig3 = q + 2 * p * cos(phi + (2*pi/3))
    //    eig2 = 3 * q - eig1 - eig3     % since trace(A) = eig1 + eig2 + eig3
    // end

    SIMD p1 = A[1] * A[1] + A[5] * A[5] + A[2] * A[2];
    SIMD q  = load1(1. / 3) * (A[0] + A[4] + A[8]);
    SIMD p2 = load1(1. / 6) * ((A[0] - q) * (A[0] - q) + (A[4] - q) * (A[4] - q) + (A[8] - q) * (A[8] - q) + load1(2) * p1);
    SIMD p  = sqrt(p2);
    SIMD rp = load1(1.) / p;
    SIMD B0 = rp * (A[0] - q);
    SIMD B1 = rp * (A[4] - q);
    SIMD B2 = rp * (A[8] - q);
    SIMD B3 = rp * A[1];
    SIMD B4 = rp * A[5];
    SIMD B5 = rp * A[2];

    SIMD r = load1(1. / 2) * (B0 * B1 * B2 + load1(2) * B3 * B4 * B5 - B5 * B5 * B1 - B3 * B3 * B2 - B4 * B4 * B0);
    SIMD phi = load1(1. / 3) * acos(r);

    eVal[0] = q + load1(2) * p * cos(phi);
    eVal[2] = q + load1(2) * p * cos(phi + load1(2) * load1(M_PI) * load1(1. / 3));
    eVal[1] = load1(3) * q - eVal[0] - eVal[2];
}

void eigSym(const SIMD A[9], SIMD eVal[3], SIMD eVec[9])
{
    eigSym(A, eVal);

    auto getVec = [] (SIMD* eVec, const SIMD A[9], const SIMD B[9]) {
        // sum up all columns to mitigate rounding errors
        SIMD v[3];
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    v[i] = v[i] + A[i * 3 + k] * B[k * 3 + j];
                }
            }
        }
        SIMD norm = rsqrt14(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        eVec[0] = v[0] * norm; eVec[1] = v[1] * norm; eVec[2] = v[2] * norm;
    };

    SIMD A0[9] = {
            A[0] - eVal[0], A[1]          , A[2]          ,
            A[1]          , A[4] - eVal[0], A[5]          ,
            A[2]          , A[5]          , A[8] - eVal[0],
    };
    SIMD A1[9] = {
            A[0] - eVal[1], A[1]          , A[2]          ,
            A[1]          , A[4] - eVal[1], A[5]          ,
            A[2]          , A[5]          , A[8] - eVal[1],
    };
    SIMD A2[9] = {
            A[0] - eVal[2], A[1]          , A[2]          ,
            A[1]          , A[4] - eVal[2], A[5]          ,
            A[2]          , A[5]          , A[8] - eVal[2],
    };
    getVec(eVec + 0, A1, A2);
    getVec(eVec + 3, A0, A2);
    getVec(eVec + 6, A0, A1);
}

void print(size_t rows, size_t cols, const SIMD *A)
{
    for (size_t s = 0; s < SIMD::size; ++s) {
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                printf("%+.4e ", A[r * cols + c][s]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

}


