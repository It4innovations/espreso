
#include "math.h"
#include "math/math.h"

#include <cstdio>

namespace espreso {

void eigSym22(const SIMD A[4], SIMD eVal[2])
{

}

void eigSym22(const SIMD A[4], SIMD eVal[2], SIMD eVec[4])
{

}

void eigSym33Desc(const SIMD A[9], SIMD eVal[3])
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

//    SIMD p1 = A[1] * A[1] + A[5] * A[5] + A[2] * A[2];
//    SIMD q  = load1(1. / 3) * (A[0] + A[4] + A[8]);
//    SIMD p2 = load1(1. / 6) * ((A[0] - q) * (A[0] - q) + (A[4] - q) * (A[4] - q) + (A[8] - q) * (A[8] - q) + load1(2) * p1);
//    SIMD p  = sqrt(p2);
//    SIMD rp = load1(1.) / p;
//    SIMD B0 = rp * (A[0] - q);
//    SIMD B1 = rp * (A[4] - q);
//    SIMD B2 = rp * (A[8] - q);
//    SIMD B3 = rp * A[1];
//    SIMD B4 = rp * A[5];
//    SIMD B5 = rp * A[2];
//
//    SIMD r = load1(1. / 2) * (B0 * B1 * B2 + load1(2) * B3 * B4 * B5 - B5 * B5 * B1 - B3 * B3 * B2 - B4 * B4 * B0);
//    SIMD phi = load1(1. / 3) * acos(r);
//
//    eVal[0] = q + load1(2) * p * cos(phi);
//    eVal[2] = q + load1(2) * p * cos(phi + load1(2) * load1(M_PI) * load1(1. / 3));
//    eVal[1] = load1(3) * q - eVal[0] - eVal[2];

    for (size_t s = 0; s < SIMD::size && !std::isnan(A[0][s]); ++s) {
        double _VAL[3], _A[9] = {
                A[0][s], A[1][s], A[2][s],
                A[3][s], A[4][s], A[5][s],
                A[6][s], A[7][s], A[8][s],
        };
        Matrix_Dense<double> _in; _in.nrows = 3; _in.ncols = 3; _in.nnz = 9; _in.vals = _A;
        Vector_Dense<double> val; val.size  = 3; val.vals = _VAL;
        math::lapack::get_eig_sym(_in, val);
        eVal[0][s] = _VAL[2]; eVal[1][s] = _VAL[1]; eVal[2][s] = _VAL[0];
    }
}

void eigSym33Asc(const SIMD A[9], SIMD eVal[3], SIMD eVec[9])
{
    for (size_t s = 0; s < SIMD::size && !std::isnan(A[0][s]); ++s) {
        double _VAL[3], _VEC[9], _A[9] = {
                A[0][s], A[1][s], A[2][s],
                A[3][s], A[4][s], A[5][s],
                A[6][s], A[7][s], A[8][s],
        };
        Matrix_Dense<double> _in; _in.nrows = 3; _in.ncols = 3; _in.nnz = 9; _in.vals = _A;
        Vector_Dense<double> val; val.size  = 3; val.vals = _VAL;
        Matrix_Dense<double> vec; vec.nrows = 3; vec.ncols = 3; vec.nnz = 9; vec.vals = _VEC;
        math::lapack::get_eig_sym(_in, val, vec);
        eVal[0][s] = _VAL[0]; eVal[1][s] = _VAL[1]; eVal[2][s] = _VAL[2];

        eVec[0][s] = _VEC[0]; eVec[1][s] = _VEC[3]; eVec[2][s] = _VEC[6];
        eVec[3][s] = _VEC[1]; eVec[4][s] = _VEC[4]; eVec[5][s] = _VEC[7];
        eVec[6][s] = _VEC[2]; eVec[7][s] = _VEC[5]; eVec[8][s] = _VEC[8];
    }
}

void eigSym33Desc(const SIMD A[9], SIMD eVal[3], SIMD eVec[9])
{
    for (size_t s = 0; s < SIMD::size && !std::isnan(A[0][s]); ++s) {
        double _VAL[3], _VEC[9], _A[9] = {
                A[0][s], A[1][s], A[2][s],
                A[3][s], A[4][s], A[5][s],
                A[6][s], A[7][s], A[8][s],
        };
        Matrix_Dense<double> _in; _in.nrows = 3; _in.ncols = 3; _in.nnz = 9; _in.vals = _A;
        Vector_Dense<double> val; val.size  = 3; val.vals = _VAL;
        Matrix_Dense<double> vec; vec.nrows = 3; vec.ncols = 3; vec.nnz = 9; vec.vals = _VEC;
        math::lapack::get_eig_sym(_in, val, vec);
        eVal[0][s] = _VAL[2]; eVal[1][s] = _VAL[1]; eVal[2][s] = _VAL[0];

        eVec[0][s] = _VEC[2]; eVec[1][s] = _VEC[5]; eVec[2][s] = _VEC[8];
        eVec[3][s] = _VEC[1]; eVec[4][s] = _VEC[4]; eVec[5][s] = _VEC[7];
        eVec[6][s] = _VEC[0]; eVec[7][s] = _VEC[3]; eVec[8][s] = _VEC[6];
    }

//    eigSym33(A, eVal);
//
//    auto getVec = [] (SIMD* eVec, const SIMD A[9], const SIMD B[9]) {
//        // sum up all columns to mitigate rounding errors
//        SIMD v[9], sum[3];
//        for (size_t i = 0; i < 3; ++i) {
//            for (size_t j = 0; j < 3; ++j) {
//                for (size_t k = 0; k < 3; ++k) {
//                    v[i * 3 + j] = v[i * 3 + j] + A[i * 3 + k] * B[k * 3 + j];
//                }
//                sum[j] = sum[j] + v[i * 3 + j] * v[i * 3 + j];
//            }
//        }
//
//        for (size_t s = 0; s < SIMD::size; ++s) {
//            for (int i = 0; i < 3; ++i) {
//                if (sum[0][s] > sum[1][s]) {
//                    if (sum[0][s] > sum[2][s]) {
//                        eVec[0][s] = v[0][s]; eVec[1][s] = v[3][s];  eVec[2][s] = v[6][s];
//                    } else {
//                        eVec[0][s] = v[2][s]; eVec[1][s] = v[5][s];  eVec[2][s] = v[8][s];
//                    }
//                } else {
//                    if (sum[1][s] > sum[2][s]) {
//                        eVec[0][s] = v[1][s]; eVec[1][s] = v[4][s];  eVec[2][s] = v[7][s];
//                    } else {
//                        eVec[0][s] = v[2][s]; eVec[1][s] = v[5][s];  eVec[2][s] = v[8][s];
//                    }
//                }
//            }
//        }
//        print(3, 3, v);
//        printf("SUM\n"); print(3, 1, sum);
//        printf("VEC\n"); print(3, 1, eVec);
//
//        SIMD norm = rsqrt14(eVec[0] * eVec[0] + eVec[1] * eVec[1] + eVec[2] * eVec[2]);
//        eVec[0] = eVec[0] * norm; eVec[1] = eVec[1] * norm; eVec[2] = eVec[2] * norm;
//    };
//
////    auto mult = [] (SIMD X[9], const SIMD A[9], const SIMD B[9]) {
////        // sum up all columns to mitigate rounding errors
////
////        for (size_t i = 0; i < 3; ++i) {
////            for (size_t j = 0; j < 3; ++j) {
////                X[i * 3 + j] = zeros();
////                for (size_t k = 0; k < 3; ++k) {
////                    X[i * 3 + j] = X[i * 3 + j] + A[i * 3 + k] * B[k * 3 + j];
////                }
////            }
////        }
////    };
//
//    SIMD A0[9] = {
//            A[0] - eVal[0], A[1]          , A[2]          ,
//            A[1]          , A[4] - eVal[0], A[5]          ,
//            A[2]          , A[5]          , A[8] - eVal[0],
//    };
//    SIMD A1[9] = {
//            A[0] - eVal[1], A[1]          , A[2]          ,
//            A[1]          , A[4] - eVal[1], A[5]          ,
//            A[2]          , A[5]          , A[8] - eVal[1],
//    };
//    SIMD A2[9] = {
//            A[0] - eVal[2], A[1]          , A[2]          ,
//            A[1]          , A[4] - eVal[2], A[5]          ,
//            A[2]          , A[5]          , A[8] - eVal[2],
//    };
//
//    getVec(eVec + 0, A1, A2);
//    getVec(eVec + 3, A0, A2);
//    getVec(eVec + 6, A0, A1);
//
//    printf("%+.12e %+.12e %+.12e\n", eVal[0][0], eVal[1][0], eVal[2][0]);
//
////    SIMD A01[9], A02[9], A12[9];
////    printf("A0\n"); print(3, 3, A0);
////    printf("A1\n"); print(3, 3, A1);
////    printf("A2\n"); print(3, 3, A2);
////
////    mult(A12, A1, A2);
////    mult(A02, A0, A2);
////    mult(A01, A0, A1);
////
////    printf("A1xA2\n"); print(3, 3, A12);
////    printf("A0xA2\n"); print(3, 3, A02);
////    printf("A0xA1\n"); print(3, 3, A01);
////
////    SIMD X[9];
////
////    mult(X, A0, A12); print(3, 3, X);
////    mult(X, A01, A2); print(3, 3, X);
//
//    printf("%+.12e %+.12e %+.12e\n", eVec[0][0], eVec[1][0], eVec[2][0]);
//    printf("%+.12e %+.12e %+.12e\n", eVec[3][0], eVec[4][0], eVec[5][0]);
//    printf("%+.12e %+.12e %+.12e\n", eVec[6][0], eVec[7][0], eVec[8][0]);
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


